#!/usr/bin/env python3
import sys, os, re
from collections import defaultdict

BASE = "/center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed"
ANNOTS = f"{BASE}/Annotations"

def detect_format(path):
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i>200: break
            if line.startswith("##gff-version") or "ID=" in line or "Parent=" in line:
                return "GFF3"
    return "GTF"

def parse_attrs_gtf(s):
    out={}
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', s):
        out[m.group(1)] = m.group(2)
    return out

def parse_attrs_gff3(s):
    out={}
    for part in s.strip().split(";"):
        if not part: continue
        if "=" in part:
            k,v = part.split("=",1)
            out[k]=v
    return out

def main(sp):
    gtf = f"{ANNOTS}/{sp}_annotation.gtf"
    fmt = detect_format(gtf)
    parse_attrs = parse_attrs_gff3 if fmt=="GFF3" else parse_attrs_gtf
    sys.stderr.write(f"[{sp}] Detected {fmt}\n")

    # Collect CDS blocks per transcript
    cds_by_tid = defaultdict(list)      # tid -> list of (chrom, start, end)
    tid_meta   = {}                     # tid -> (gid, gname, chrom, strand)

    with open(gtf) as fh:
        for line in fh:
            if not line or line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f)<9: continue
            chrom, source, feature, start, end, score, strand, frame, attr = f
            start = int(start); end = int(end)
            A = parse_attrs(attr)

            if feature in ("transcript","mRNA","gene"):
                # record metadata when available
                if fmt=="GFF3":
                    if feature in ("mRNA","transcript"):
                        tid = A.get("ID")
                        gid = A.get("Parent") or A.get("gene_id")
                    elif feature=="gene":
                        tid = None
                        gid = A.get("ID")
                    gname = A.get("gene_name") or A.get("Name") or A.get("gene")
                else:
                    tid = A.get("transcript_id")
                    gid = A.get("gene_id")
                    gname = A.get("gene_name") or A.get("gene_symbol") or A.get("Name")
                if tid:
                    tid_meta.setdefault(tid, (gid, gname, chrom, strand))
                continue

            if feature!="CDS":
                continue

            if fmt=="GFF3":
                parent = A.get("Parent", "")
                tid = parent.split(",")[0] if parent else A.get("transcript_id")
                gid = A.get("gene_id") or A.get("gene")
                gname = A.get("gene_name") or A.get("Name") or A.get("gene")
            else:
                tid = A.get("transcript_id")
                gid = A.get("gene_id")
                gname = A.get("gene_name")

            if not tid:  # synthesize if only gene present
                if gid:
                    tid = f"{gid}.tx"
                else:
                    continue

            # ensure meta
            if tid not in tid_meta:
                tid_meta[tid] = (gid, gname, chrom, strand if strand in ("+","-") else "+")
            cds_by_tid[tid].append((chrom, start, end))

    # Write BED12 (one line per transcript; CDS blocks only)
    bed12 = f"{sp}_cds.bed12"
    idx   = f"{sp}_cds_index.tsv"
    with open(bed12, "w") as fb, open(idx, "w") as fi:
        fi.write("gene_id\ttranscript_id\tgene_name\tstrand\texon_count\tcds_len\n")
        for tid, blocks in cds_by_tid.items():
            gid, gname, chrom, strand = tid_meta.get(tid, (None, "", blocks[0][0], "+"))
            blocks.sort(key=lambda x: x[1])  # ascending by start
            chromStart = blocks[0][1]-1      # 0-based
            chromEnd   = blocks[-1][2]
            blockSizes = []
            blockStarts= []
            cds_len = 0
            for (c,s,e) in blocks:
                sz = e - s + 1
                cds_len += sz
                blockSizes.append(str(sz))
                blockStarts.append(str(s-1 - chromStart))
            name = f"{tid}|gene:{gid}|gene_name:{'' if gname is None else gname}"
            score = 0
            thickStart = chromStart
            thickEnd   = chromEnd
            itemRgb = 0
            blockCount = len(blockSizes)
            fb.write("\t".join(map(str,[
                chrom, chromStart, chromEnd, name, score, strand,
                thickStart, thickEnd, itemRgb, blockCount,
                ",".join(blockSizes)+",",
                ",".join(blockStarts)+","
            ]))+"\n")
            fi.write("\t".join([str(gid or ""), str(tid), str(gname or ""), strand, str(blockCount), str(cds_len)])+"\n")

    sys.stderr.write(f"[{sp}] Wrote {bed12} and {idx}\n")

if __name__=="__main__":
    if len(sys.argv)!=2:
        sys.exit("Usage: gtf_to_bed12_cds.py <SPECIES_CODE>")
    main(sys.argv[1])

