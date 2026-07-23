#!/usr/bin/env python3
import sys, os, subprocess, re
from collections import defaultdict

BASE = "/center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed"
GENOMES = f"{BASE}/Genomes"
ANNOTS  = f"{BASE}/Annotations"

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
def rc(s): return s.translate(COMP)[::-1]

def faidx_seq(genome_fa, chrom, start, end):
    region = f"{chrom}:{start}-{end}"
    try:
        r = subprocess.run(["samtools","faidx", genome_fa, region],
                           capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError:
        return ""  # contig mismatch or bad coords
    return "".join([l.strip() for l in r.stdout.splitlines() if not l.startswith(">")])

def guess_format_and_attr_parser(gtf_path):
    # Peek 200 lines to decide GTF vs GFF3 and build attribute parser
    is_gff3 = False
    with open(gtf_path) as fh:
        for i, line in enumerate(fh):
            if i > 200: break
            if line.startswith("##gff-version") or "ID=" in line or "Parent=" in line:
                is_gff3 = True
                break
    if not is_gff3:
        # GTF-style: key "gene_id" "transcript_id" "gene_name"
        def parse_attrs(s):
            out={}
            for m in re.finditer(r'(\S+)\s+"([^"]+)"', s):
                out[m.group(1)] = m.group(2)
            return out
        return "GTF", parse_attrs
    else:
        # GFF3-style: semicolons, key=value
        def parse_attrs(s):
            out={}
            for part in s.strip().split(";"):
                if not part: continue
                if "=" in part:
                    k,v = part.split("=",1)
                    out[k]=v
            return out
        return "GFF3", parse_attrs

def main(sp):
    genome = f"{GENOMES}/{sp}_genome.fa"
    gtf    = f"{ANNOTS}/{sp}_annotation.gtf"
    if not os.path.exists(genome): sys.exit(f"Missing genome: {genome}")
    if not os.path.exists(gtf):    sys.exit(f"Missing annotation: {gtf}")
    # index genome if needed
    if not os.path.exists(genome + ".fai"):
        subprocess.run(["samtools","faidx", genome], check=True)

    fmt, parse_attrs = guess_format_and_attr_parser(gtf)
    print(f"[{sp}] Detected {fmt} attributes", file=sys.stderr)

    # collect contigs from .fai to warn about mismatches
    fai_contigs = set()
    with open(genome + ".fai") as fh:
        for line in fh:
            fai_contigs.add(line.split("\t")[0])

    # transcript metadata and CDS blocks
    tx_meta = {}  # tid -> {gid, gname, chrom, strand}
    cds_blocks = defaultdict(list)  # tid -> [(start,end)]
    cds_count = 0

    with open(gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attr = f
            start = int(start); end = int(end)
            A = parse_attrs(attr)

            # map IDs depending on format
            if fmt == "GTF":
                gid = A.get("gene_id")
                tid = A.get("transcript_id")
                gnm = A.get("gene_name") or A.get("gene_symbol") or A.get("Name")
            else:
                # GFF3: gene has ID=geneX, transcript has ID=txY Parent=geneX; CDS has Parent=txY (often)
                gid = A.get("gene_id")  # rare in GFF3
                if feature == "gene":
                    gid = A.get("ID", gid)
                if feature in ("mRNA","transcript"):
                    tid = A.get("ID")
                    if not gid:
                        gid = A.get("Parent") or A.get("gene")
                elif feature == "CDS":
                    # CDS Parent should be transcript ID; sometimes comma-separated
                    parent = A.get("Parent", "")
                    tid = parent.split(",")[0] if parent else None
                else:
                    tid = A.get("transcript_id") or A.get("ID")

                # try to recover gene name
                gnm = A.get("gene_name") or A.get("Name") or A.get("gene")

            if feature in ("transcript","mRNA","gene"):
                if tid:
                    d = tx_meta.get(tid, {"gene_id":gid, "gene_name":gnm, "chrom":chrom, "strand":strand})
                    if gid: d["gene_id"]=gid
                    if gnm: d["gene_name"]=gnm
                    d["chrom"]=chrom
                    d["strand"]=strand if strand in ("+","-") else d.get("strand","+")
                    tx_meta[tid]=d
                continue

            if feature != "CDS":
                continue
            # CDS features
            if fmt == "GTF":
                tid = A.get("transcript_id")
                gid = A.get("gene_id") or gid
                gnm = A.get("gene_name") or gnm
            # ensure we have a transcript id
            if not tid:
                # last resort: some GFF3s put CDS Parent=gene; then we synthesize a per-gene transcript
                if gid:
                    tid = f"{gid}.tx"
                else:
                    continue
            # fill transcript meta if missing
            if tid not in tx_meta:
                tx_meta[tid] = {"gene_id": gid, "gene_name": gnm, "chrom": chrom, "strand": strand if strand in ("+","-") else "+"}
            cds_blocks[tid].append((chrom, start, end))
            cds_count += 1

    print(f"[{sp}] Parsed CDS features: {cds_count}", file=sys.stderr)
    if cds_count == 0:
        print(f"[{sp}] WARNING: no CDS features found. Check file format & attributes.", file=sys.stderr)

    out_tx = f"{sp}_cds_transcripts.fa"
    out_gene = f"{sp}_cds_longest_per_gene.fa"
    out_idx = f"{sp}_cds_index.tsv"

    gene_best = {}  # gid -> (len, tid, seq, gname, exon_count)
    wrote = 0

    with open(out_tx, "w") as ftx, open(out_idx, "w") as findx:
        findx.write("gene_id\ttranscript_id\tgene_name\tstrand\texon_count\tcds_len\n")
        for i, (tid, blocks) in enumerate(cds_blocks.items(), 1):
            meta = tx_meta.get(tid, {})
            gid  = meta.get("gene_id") or tid
            gnm  = meta.get("gene_name") or ""
            strand = meta.get("strand","+")
            # order blocks by genomic coord; reverse later if minus strand
            blocks.sort(key=lambda x: x[1])
            seq_parts = []
            bad_contig = False
            for chrom, s, e in blocks:
                if chrom not in fai_contigs:
                    bad_contig = True
                    break
                seq_parts.append(faidx_seq(genome, chrom, s, e))
            if bad_contig or any(len(s)==0 for s in seq_parts):
                # contig mismatch or empty fetch; skip this transcript
                continue
            seq = "".join(seq_parts)
            if strand == "-":
                seq = rc(seq)
            if not seq:
                continue
            # write transcript
            ftx.write(f">{tid}|gene:{gid}|gene_name:{gnm}\n")
            for j in range(0,len(seq),60):
                ftx.write(seq[j:j+60]+"\n")
            exon_count = len(blocks)
            findx.write(f"{gid}\t{tid}\t{gnm}\t{strand}\t{exon_count}\t{len(seq)}\n")
            wrote += 1
            # longest per gene
            prev = gene_best.get(gid)
            cur  = (len(seq), tid, seq, gnm, exon_count, strand)
            if (prev is None) or (len(seq) > prev[0]):
                gene_best[gid] = cur
            if i % 1000 == 0:
                print(f"[{sp}] processed transcripts: {i} (written {wrote})", file=sys.stderr)

    with open(out_gene, "w") as fg:
        for gid, (L, tid, seq, gnm, exon_count, strand) in gene_best.items():
            fg.write(f">{gid}|rep_tid:{tid}|gene_name:{gnm}|len:{L}\n")
            for j in range(0,len(seq),60):
                fg.write(seq[j:j+60]+"\n")

    print(f"[{sp}] DONE: wrote {wrote} transcripts, {len(gene_best)} genes", file=sys.stderr)
    print(f"[{sp}] Outputs: {out_tx}, {out_gene}, {out_idx}", file=sys.stderr)

if __name__ == "__main__":
    if len(sys.argv)!=2:
        sys.exit("Usage: extract_cds_robust.py <SPECIES_CODE>")
    main(sys.argv[1])

