# gib_extract_accessions.py
import csv, os, gzip, re, sys
from collections import defaultdict

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def parse_attributes(attr):
    """Parse GTF key "value" and GFF3 key=value; returns dict."""
    d = {}
    for k, v in re.findall(r'(\S+)\s+"([^"]+)"', attr):  # GTF
        d[k] = v
    for k, v in re.findall(r'(\S+)=([^;]+)', attr):      # GFF3
        d[k] = v
    return d

def first_present(attrs, keys):
    for k in keys:
        v = attrs.get(k)
        if v and v not in (".", "-"):
            return v
    return None

# Regexes for common protein accessions (RefSeq, UniProt) + generic NCBI-like (e.g., KAI4832453.1)
re_refseq = re.compile(r'\b((?:NP|XP|YP|WP)_[0-9]+(?:\.[0-9]+)?)\b')
re_uniprot = re.compile(r'\b((?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]{5})(?:-\d+)?)\b')
re_ncbi_generic = re.compile(r'\b([A-Z]{3}[A-Z0-9]{5,}(?:\.[0-9]+)?)\b')  # e.g., KAI4832453.1

def harvest_protein_ids(attrs):
    """Collect protein accessions from multiple places/keys."""
    hits = set()

    # direct keys seen in various GTF/GFFs
    for key in ("protein_id", "protein", "protein_accession", "peptide_id", "translation_id"):
        v = attrs.get(key)
        if v and v not in (".", "-"):
            # Some files have multiple IDs separated by commas/spaces
            for tok in re.split(r'[,\s]+', v):
                if tok:
                    hits.add(tok)

    # Dbxref can hold RefSeq/UniProt, e.g., Dbxref=RefSeq:XP_...,UniProtKB/TrEMBL:Q9...,Genbank:...
    dbx = attrs.get("Dbxref") or attrs.get("db_xref") or attrs.get("DbxRef")
    if dbx:
        # split on comma/semicolon
        for part in re.split(r'[;,]', dbx):
            # take text after colon if present (RefSeq:XP_..., UniProtKB/Swiss-Prot:P12345)
            if ":" in part:
                part = part.split(":", 1)[1]
            part = part.strip()
            if not part:
                continue
            # run regex matchers
            for rx in (re_refseq, re_uniprot, re_ncbi_generic):
                for m in rx.findall(part):
                    hits.add(m)

    return hits

def extract_gib_accessions(gtf_file, out_wide, out_long):
    # Pass 1: collect gene↔name, transcript↔gene, ID↔Parent (for GFF3 chains)
    gene_name_by_gene = {}         # gene key (gene_id or gene-feature ID) -> gene_name
    transcript_to_gene = {}        # transcript ID -> gene key
    id_to_parent = {}              # ID -> Parent (for climbing)
    feature_type = {}              # ID -> ftype (gene/mRNA/transcript)

    with open_text(gtf_file) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2]
            attrs = parse_attributes(f[8])

            # Track generic IDs and parents (GFF3)
            feat_id = first_present(attrs, ["ID", "transcript_id"])
            parent  = first_present(attrs, ["Parent"])
            if feat_id:
                feature_type[feat_id] = ftype
            if feat_id and parent:
                id_to_parent[feat_id] = parent

            # Gene name capture from any line that provides it
            gene_id   = first_present(attrs, ["gene_id", "gene"])    # gene key
            gene_name = first_present(attrs, ["gene_name", "Name", "gene_symbol"])
            if ftype == "gene":
                gene_feat_id = first_present(attrs, ["ID"]) or gene_id
                if gene_name and gene_feat_id:
                    gene_name_by_gene[gene_feat_id] = gene_name
                if gene_name and gene_id:
                    gene_name_by_gene[gene_id] = gene_name
            else:
                if gene_name and gene_id:
                    gene_name_by_gene[gene_id] = gene_name

            # transcript -> gene mapping
            if ftype in ("transcript", "mRNA"):
                tid = first_present(attrs, ["transcript_id", "ID"])
                gid = gene_id or parent
                if tid and gid:
                    transcript_to_gene[tid] = gid

    # Pass 2: scan ALL features for protein accessions; resolve to gene symbol via maps
    gene_to_proteins = defaultdict(set)
    total_lines = 0
    protein_lines = 0
    resolved = 0

    def resolve_symbol(attrs):
        # Direct symbol on the same line?
        sym = first_present(attrs, ["gene_name", "Name", "gene_symbol"])
        if sym:
            return sym
        # Explicit gene on this line → lookup
        gid = first_present(attrs, ["gene_id", "gene"])
        if gid and gid in gene_name_by_gene:
            return gene_name_by_gene[gid]
        # Climb via transcript/Parent chain
        tid = first_present(attrs, ["transcript_id", "Parent"])
        if tid:
            # transcript -> gene
            gid2 = transcript_to_gene.get(tid) or id_to_parent.get(tid)
            if gid2 and gid2 in gene_name_by_gene:
                return gene_name_by_gene[gid2]
            # Sometimes transcript carries a gene_name mapping we captured
            if tid in gene_name_by_gene:
                return gene_name_by_gene[tid]
        # As last resort, keep a stable gene-like ID so pairs aren’t lost
        return gid or tid or "UNRESOLVED"

    with open_text(gtf_file) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            total_lines += 1
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = parse_attributes(f[8])

            prots = harvest_protein_ids(attrs)
            if not prots:
                continue
            protein_lines += 1

            symbol = resolve_symbol(attrs)
            for p in prots:
                gene_to_proteins[symbol].add(p)
                if symbol != "UNRESOLVED":
                    resolved += 1

    # Write outputs
    os.makedirs(os.path.dirname(out_wide) or ".", exist_ok=True)

    with open(out_wide, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["Gene_Symbol", "Protein_Accession_IDs"])
        for gene, prots in sorted(gene_to_proteins.items()):
            w.writerow([gene, ",".join(sorted(prots))])

    with open(out_long, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["Gene_Symbol", "Protein_Accession_ID"])
        for gene, prots in sorted(gene_to_proteins.items()):
            for p in sorted(prots):
                w.writerow([gene, p])

    # Progress
    print(f"[INFO] total lines: {total_lines}", file=sys.stderr)
    print(f"[INFO] lines with protein accessions found: {protein_lines}", file=sys.stderr)
    print(f"[INFO] mapped pairs (non-UNRESOLVED): {resolved}", file=sys.stderr)
    if not gene_to_proteins:
        print("[WARN] No protein accessions were discovered. "
              "GIB GTF may lack Dbxref/protein IDs; consider using gffread to make CDS FASTA and skip accessions.",
              file=sys.stderr)

if __name__ == "__main__":
    gtf_in = "GIB_annotation.gtf"        # or .gtf.gz
    out_dir = "./protein_mappings"
    os.makedirs(out_dir, exist_ok=True)
    out_wide = f"{out_dir}/GIB_protein_mapping.tsv"
    out_long = f"{out_dir}/GIB_protein_pairs.tsv"
    print(f"Processing GIB from {gtf_in} ...")
    extract_gib_accessions(gtf_in, out_wide, out_long)
    print(f"Saved:\n  {out_wide}\n  {out_long}")
