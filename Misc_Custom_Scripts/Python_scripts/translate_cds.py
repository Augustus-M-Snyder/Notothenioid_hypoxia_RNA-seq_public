#!/usr/bin/env python3
import sys, argparse
from textwrap import wrap

# Standard (NCBI table 1) translate dict; ambiguous codons -> 'X'
CODON = {
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
 'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
 'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
 'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
 'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def translate(seq):
    seq = seq.upper().replace('U','T')
    prot = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) < 3: break
        aa = CODON.get(codon, 'X')
        prot.append(aa)
    p = ''.join(prot)
    # trim single trailing stop; keep internal ones for QC
    return p[:-1] if p.endswith('*') else p

def read_fa(path):
    with open(path) as fh:
        hdr, seq = None, []
        for line in fh:
            if line.startswith('>'):
                if hdr is not None:
                    yield hdr, ''.join(seq)
                hdr, seq = line.strip()[1:], []
            else:
                seq.append(line.strip())
        if hdr is not None:
            yield hdr, ''.join(seq)

def write_fa(path, records):
    with open(path, 'w') as out:
        for h, s in records:
            out.write(f'>{h}\n')
            for chunk in wrap(s, 60):
                out.write(chunk + '\n')

def main():
    ap = argparse.ArgumentParser(description="Translate CDS FASTA to protein FASTA (QC included).")
    ap.add_argument("cds_fa", help="*_cds_transcripts.fa or *_cds_longest_per_gene.fa")
    ap.add_argument("-o","--out", help="output protein FASTA (default: auto)")
    ap.add_argument("--qc", help="QC report TSV (default: auto)")
    args = ap.parse_args()

    out = args.out or args.cds_fa.replace(".fa", ".prot.fa")
    qc  = args.qc  or args.cds_fa.replace(".fa", ".translate.qc.tsv")

    kept, drops = [], []
    internal_stop = 0
    total = 0

    for hdr, nt in read_fa(args.cds_fa):
        total += 1
        if len(nt) % 3 != 0:
            drops.append((hdr, len(nt), "len%3!=0"))
            continue
        p = translate(nt)
        if '*' in p:
            internal_stop += 1
        kept.append((hdr, p))

    write_fa(out, kept)
    with open(qc, "w") as fh:
        fh.write("id\tnt_len\tflag\n")
        for h,l,flag in drops:
            fh.write(f"{h}\t{l}\t{flag}\n")
    sys.stderr.write(
        f"[translate] input={args.cds_fa} out={out} kept={len(kept)} "
        f"dropped={len(drops)} internal_stop_proteins={internal_stop}\n"
    )

if __name__ == "__main__":
    main()

