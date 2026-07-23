#!/usr/bin/env python3
import sys, os, subprocess, tempfile, shutil
from collections import defaultdict, namedtuple

BASE = "/center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed"
GENOMES = f"{BASE}/Genomes"
ANNOTS  = f"{BASE}/Annotations"

Comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def rc(seq): return seq.translate(Comp)[::-1]

def faidx_seq(genome_fa, chrom, start, end, strand):
    # 1-based inclusive coords expected for faidx
    region = f"{chrom}:{start}-{end}"
    r = subprocess.run(["samtools","faidx", genome_fa, region],
                       capture_output=True, text=True, check=True)
    seq = "".join([l.strip() for l in r.stdout.splitlines() if not l.startswith(">")])
    if strand == "-":
        seq = rc(seq)
    return seq

def parse_gtf_collect_cds(gtf_path):
    # returns: dict[transcript_id] -> dict with gene_id, gene_name, chrom, strand, list of (start,end) CDS exons
    Tx = {}
    CDS = defaultdict(list)
    # Basic attribute parser
    def attrs2dict(s):
        out={}
        for part in s.strip().split(";"):
            part=part.strip()
            if not part: continue
            if " " in part:
                k,v = part.split(" ",1)
                v = v.strip().strip('"')
                out[k]=v
        return out

    with open(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9: continue
            chrom, source, feature, start, end, score, strand, frame, attr = f
            if feature != "CDS": 
                # keep gene/transcript metadata when we see them
                if feature in ("transcript","mRNA","exon","gene"):
                    d = attrs2dict(attr)
                    gid = d.get("gene_id")
                    tid = d.get("transcript_id")
                    gnm = d.get("gene_name") or d.get("gene_symbol") or d.get("Name")
                    if tid:
                        info = Tx.get(tid, {"gene_id":gid,"gene_name":gnm,"chrom":chrom,"strand":strand})
                        if gid: info["gene_id"]=gid
                        if gnm: info["gene_name"]=gnm
                        info["chrom"]=chrom
                        info["strand"]=strand if strand in ("+","-") else info.get("strand","+")
                        Tx[tid]=info
                continue
            # CDS rows
            d = attrs2dict(attr)
            tid = d.get("transcript_id")
            gid = d.get("gene_id")
            gnm = d.get("gene_name") or d.get("gene_symbol") or d.get("Name")
            strand = strand if strand in ("+","-") else "+"
            if not tid: 
                continue
            if tid not in Tx:
                Tx[tid] = {"gene_id":gid,"gene_name":gnm,"chrom":chrom,"strand":strand}
            else:
                if gid: Tx[tid]["gene_id"]=gid
                if gnm: Tx[tid]["gene_name"]=gnm
                Tx[tid]["chrom"]=chrom
                Tx[tid]["strand"]=strand
            s = int(start); e = int(end)
            CDS[tid].append((s,e))
    # sort exons along transcription orientation
    for tid in CDS:
        strand = Tx.get(tid,{}).get("strand","+")
        CDS[tid].sort(key=lambda x: x[0], reverse=(strand=="-"))
    return Tx, CDS

def main(sp):
    genome = f"{GENOMES}/{sp}_genome.fa"
    gtf    = f"{ANNOTS}/{sp}_annotation.gtf"
    assert os.path.exists(genome), f"Missing genome: {genome}"
    assert os.path.exists(gtf),    f"Missing GTF: {gtf}"

    # Ensure faidx
    if not os.path.exists(genome + ".fai"):
        subprocess.run(["samtools","faidx", genome], check=True)

    Tx, CDS = parse_gtf_collect_cds(gtf)

    out_tx_fa = f"{sp}_cds_transcripts.fa"
    out_gene_fa = f"{sp}_cds_longest_per_gene.fa"
    out_idx = f"{sp}_cds_index.tsv"

    # Gather sequences per transcript
    gene2best = {}  # gene_id -> (cds_len, tid, seq, gene_name, strand, exon_count)
    with open(out_tx_fa,"w") as ft, open(out_idx,"w") as fi:
        fi.write("\t".join(["gene_id","transcript_id","gene_name","strand","exon_count","cds_len"])+"\n")
        for tid, blocks in CDS.items():
            meta = Tx.get(tid, {})
            gid  = meta.get("gene_id") or tid  # fallback
            gnm  = meta.get("gene_name") or ""
            chrom= meta.get("chrom")
            strand= meta.get("strand","+")
            if not chrom:
                continue
            seq_parts=[]
            for s,e in blocks:
                seq_parts.append(faidx_seq(genome, chrom, s, e, strand))
            cds_seq = "".join(seq_parts)
            if not cds_seq: 
                continue
            hdr = f">{tid}|gene:{gid}|gene_name:{gnm}"
            ft.write(hdr+"\n")
            # wrap 60
            for i in range(0,len(cds_seq),60):
                ft.write(cds_seq[i:i+60]+"\n")
            exon_count = len(blocks)
            cds_len = len(cds_seq)
            fi.write("\t".join([gid, tid, gnm, strand, str(exon_count), str(cds_len)])+"\n")
            # longest per gene
            prev = gene2best.get(gid)
            cur  = (cds_len, tid, cds_seq, gnm, strand, exon_count)
            if (prev is None) or (cds_len > prev[0]):
                gene2best[gid] = cur

    # Write longest CDS per gene
    with open(out_gene_fa,"w") as fg:
        for gid,(cds_len, tid, seq, gnm, strand, exon_count) in gene2best.items():
            hdr = f">{gid}|rep_tid:{tid}|gene_name:{gnm}|len:{cds_len}"
            fg.write(hdr+"\n")
            for i in range(0,len(seq),60):
                fg.write(seq[i:i+60]+"\n")

    print(f"[OK] {sp}: {out_tx_fa}, {out_gene_fa}, {out_idx}")

if __name__ == "__main__":
    if len(sys.argv)!=2:
        print("Usage: extract_cds_no_gffread.py <SPECIES_CODE>", file=sys.stderr); sys.exit(1)
    main(sys.argv[1])

