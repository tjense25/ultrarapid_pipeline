#!/usr/bin/env python3

import argparse, csv, gzip, re, subprocess, sys, tempfile
from collections import defaultdict
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument("--cnv", required=True)
p.add_argument("--clinvar", required=True)
p.add_argument("--dgv", required=True)
p.add_argument("--collins", required=True)
p.add_argument("--mane", required=True)
p.add_argument("--genes", required=True)
p.add_argument("--output", required=True)
p.add_argument("--min-size", type=int, default=50000)
p.add_argument("--min-probes", type=int, default=10)
p.add_argument("--reciprocal-overlap", type=float, default=0.5)
p.add_argument("--clinvar-reference-fraction", type=float, default=0.9)
p.add_argument("--dgv-patient-fraction", type=float, default=0.9)
p.add_argument("--dgv-common-af", type=float, default=0.001)
p.add_argument("--phaplo", type=float, default=0.86)
p.add_argument("--ptriplo", type=float, default=0.94)
p.add_argument("--min-cds-fraction", type=float, default=0.84)
args = p.parse_args()

def op(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def clean(value):
    value = str(value or "").replace("\t", " ").replace("\n", " ").strip()
    return value or "."

def svtype(value):
    value = value.lower()
    if "loss" in value or value in {"del", "deletion"}:
        return "LOSS"
    if "gain" in value or value in {"dup", "duplication"}:
        return "GAIN"
    return None

def overlap_fractions(a1, a2, b1, b2):
    overlap = max(0, min(a2, b2) - max(a1, b1))
    if not overlap:
        return 0, 0, 0
    patient = overlap / (a2 - a1)
    reference = overlap / (b2 - b1)
    return patient, reference, min(patient, reference)

def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return merged

def intersect(a, b, *options):
    result = subprocess.run(
        ["bedtools", "intersect", "-a", str(a), "-b", str(b), *options],
        text=True,
        capture_output=True,
    )
    if result.returncode:
        raise SystemExit(result.stderr)
    return result.stdout.splitlines()

def chrom_key(chrom):
    c = chrom.removeprefix("chr")
    return (0, int(c)) if c.isdigit() else (1, {"X": 23, "Y": 24}.get(c, 100), c)

disease_genes = set()
with op(args.genes) as f:
    for line in f:
        gene = line.strip().split()[0] if line.strip() else ""
        if gene and not gene.startswith("#"):
            disease_genes.add(gene)

scores = {}
with op(args.collins) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        x = line.rstrip().split("\t")
        scores[x[0]] = float(x[1]), float(x[2])

cnvs = {}
with op(args.cnv) as f:
    header = next(line for line in f if line.strip()).rstrip().lstrip("#").split("\t")
    reader = csv.DictReader(f, fieldnames=header, delimiter="\t")

    for row in reader:
        row = {k.lower(): v for k, v in row.items()}
        chrom = row["chromosome"]
        start, end = int(row["start"]), int(row["end"])
        cn = float(row["cn"])
        probes = float(row["probes"])

        if cn == 2 or probes < args.min_probes or end - start < args.min_size or chrom.removeprefix("chr") in {"X", "Y"}:
            continue

        qid = f"CNV{len(cnvs) + 1}"
        cnvs[qid] = {
            "chrom": chrom,
            "start": start,
            "end": end,
            "cn": cn,
            "log2": row["log2"],
            "depth": row["depth"],
            "type": "LOSS" if cn < 2 else "GAIN",
        }

if not cnvs:
    raise SystemExit("No CNVs remained after basic filtering")

query_chroms = {x["chrom"] for x in cnvs.values()}

cds_raw = defaultdict(list)
with op(args.mane) as f:
    for line in f:
        if line.startswith("#"):
            continue
        x = line.rstrip().split("\t")

        if len(x) != 9 or x[2] != "CDS" or x[0] not in query_chroms or 'tag "MANE_Select"' not in x[8]:
            continue

        attrs = dict(re.findall(r'(\S+) "([^"]+)"', x[8]))
        gene = attrs.get("gene_name")

        if gene in disease_genes:
            cds_raw[(x[0], gene)].append((int(x[3]) - 1, int(x[4])))

cds_meta, cds_intervals = {}, {}
for i, (key, intervals) in enumerate(cds_raw.items(), 1):
    token = f"G{i}"
    cds_meta[token] = key
    cds_intervals[token] = merge_intervals(intervals)

clinvar = {}
with op(args.clinvar) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        x = line.rstrip().split("\t")

        if len(x) < 11 or x[0] not in query_chroms:
            continue

        start, end = int(x[1]), int(x[2])
        if end - start < args.min_size:
            continue

        significance = x[9].lower()
        if "pathogenic" not in significance or any(v in significance for v in ["benign", "uncertain", "conflicting"]):
            continue

        ctype = svtype(x[10])
        if not ctype:
            continue

        rid = f"C{len(clinvar) + 1}"
        clinvar[rid] = {
            "chrom": x[0],
            "start": start,
            "end": end,
            "id": clean(x[3]),
            "type": ctype,
            "phenotype": clean(x[18] if len(x) > 18 else "."),
            "review": clean(x[21] if len(x) > 21 else "."),
        }

dgv = {}
with op(args.dgv) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        x = line.rstrip().split("\t")

        if len(x) < 15 or x[0] not in query_chroms:
            continue

        dtype = svtype(x[14])
        if not dtype:
            continue

        start, end = int(x[1]), int(x[2])

        try:
            sizes = [int(v) for v in x[10].rstrip(",").split(",")]
            starts = [int(v) for v in x[11].rstrip(",").split(",")]
            i = max(range(min(len(sizes), len(starts))), key=lambda j: sizes[j])
            start, end = start + starts[i], start + starts[i] + sizes[i]
        except (ValueError, IndexError):
            pass

        af = None
        if len(x) > 26:
            try:
                af = float(x[26].replace("%", "")) / 100
            except ValueError:
                pass

        rid = f"D{len(dgv) + 1}"
        dgv[rid] = {
            "chrom": x[0],
            "start": start,
            "end": end,
            "id": clean(x[3]),
            "type": dtype,
            "af": af,
        }

clinvar_matches = defaultdict(list)
dgv_matches = defaultdict(list)
cds_overlap = defaultdict(lambda: defaultdict(int))

with tempfile.TemporaryDirectory() as tmp:
    tmp = Path(tmp)
    query_bed = tmp / "query.bed"
    clinvar_bed = tmp / "clinvar.bed"
    dgv_bed = tmp / "dgv.bed"
    cds_bed = tmp / "cds.bed"

    with open(query_bed, "w") as out:
        for qid, x in cnvs.items():
            out.write(f'{x["chrom"]}\t{x["start"]}\t{x["end"]}\t{qid}\t{x["type"]}\n')

    with open(clinvar_bed, "w") as out:
        for rid, x in clinvar.items():
            out.write(f'{x["chrom"]}\t{x["start"]}\t{x["end"]}\t{rid}\t{x["type"]}\n')

    with open(dgv_bed, "w") as out:
        for rid, x in dgv.items():
            out.write(f'{x["chrom"]}\t{x["start"]}\t{x["end"]}\t{rid}\t{x["type"]}\n')

    with open(cds_bed, "w") as out:
        for token, intervals in cds_intervals.items():
            chrom = cds_meta[token][0]
            for start, end in intervals:
                out.write(f"{chrom}\t{start}\t{end}\t{token}\n")

    if clinvar:
        for line in intersect(query_bed, clinvar_bed, "-wa", "-wb"):
            x = line.split("\t")
            qid, rid = x[3], x[8]

            if cnvs[qid]["type"] != clinvar[rid]["type"]:
                continue

            match = dict(clinvar[rid])
            patient, reference, reciprocal = overlap_fractions(
                cnvs[qid]["start"], cnvs[qid]["end"],
                match["start"], match["end"],
            )
            match.update(patient_fraction=patient, reference_fraction=reference, reciprocal=reciprocal)

            if reciprocal >= args.reciprocal_overlap or reference >= args.clinvar_reference_fraction:
                clinvar_matches[qid].append(match)

    if dgv:
        for line in intersect(query_bed, dgv_bed, "-wa", "-wb"):
            x = line.split("\t")
            qid, rid = x[3], x[8]

            if cnvs[qid]["type"] != dgv[rid]["type"]:
                continue

            match = dict(dgv[rid])
            patient, reference, reciprocal = overlap_fractions(
                cnvs[qid]["start"], cnvs[qid]["end"],
                match["start"], match["end"],
            )
            match.update(patient_fraction=patient, reference_fraction=reference, reciprocal=reciprocal)
            dgv_matches[qid].append(match)

    if cds_meta:
        for line in intersect(query_bed, cds_bed, "-wo"):
            x = line.split("\t")
            cds_overlap[x[3]][x[8]] += int(x[9])

rows = []
counts = {1: 0, 2: 0}
dgv_excluded = 0

for qid, cnv in cnvs.items():
    common_dgv = [
        x for x in dgv_matches[qid]
        if x["patient_fraction"] >= args.dgv_patient_fraction
        and x["af"] is not None
        and x["af"] > args.dgv_common_af
    ]

    if common_dgv:
        dgv_excluded += 1
        continue

    overlapping_genes = []
    for token, overlap_bp in cds_overlap[qid].items():
        chrom, gene = cds_meta[token]
        total_cds = sum(end - start for start, end in cds_intervals[token])
        fraction = overlap_bp / total_cds

        score = None
        if gene in scores:
            phaplo, ptriplo = scores[gene]
            score = phaplo if cnv["type"] == "LOSS" else ptriplo

        overlapping_genes.append({
            "gene": gene,
            "fraction": fraction,
            "score": score,
        })

    overlapping_genes.sort(
        key=lambda x: (
            x["score"] is None,
            -(x["score"] if x["score"] is not None else 0),
            x["gene"],
        )
    )

    threshold = args.phaplo if cnv["type"] == "LOSS" else args.ptriplo
    qualifying_genes = [
        x for x in overlapping_genes
        if x["fraction"] >= args.min_cds_fraction
        and x["score"] is not None
        and x["score"] >= threshold
    ]

    if clinvar_matches[qid]:
        tier, reason = 1, "CLINVAR_PLP"
        reported_genes = overlapping_genes
    elif qualifying_genes:
        tier, reason = 2, "DOSAGE_SENSITIVE_DISEASE_GENE"
        reported_genes = qualifying_genes
    else:
        continue

    best_clinvar = None
    if clinvar_matches[qid]:
        best_clinvar = max(
            clinvar_matches[qid],
            key=lambda x: (
                x["reciprocal"] >= args.reciprocal_overlap,
                x["reciprocal"] if x["reciprocal"] >= args.reciprocal_overlap else x["reference_fraction"],
                x["reference_fraction"],
            ),
        )

    best_dgv = None
    if dgv_matches[qid]:
        best_dgv = max(
            dgv_matches[qid],
            key=lambda x: (
                x["patient_fraction"],
                -1 if x["af"] is None else x["af"],
            ),
        )

    chrom=cnv["chrom"]
    original_start = int(cnv["start"])
    original_end = int(cnv["end"])
    coords = (
        f"{chrom}:"
        f"{original_start}-"
        f"{original_end}"
    )
    cnv_length = original_end - original_start
    slop = min(round(cnv_length/2),1000000)
    report_start = max(0, original_start - slop)
    report_end = original_end + slop

    length_string = "%.2fMb" % (cnv_length/1e6) if cnv_length>=1e6 else "%.1fKb" % (cnv_length/1e3)

    rows.append({
        "CHROM": cnv["chrom"],
        "START": report_start,
        "END": report_end,
        "COORDS": coords,
        "LENGTH": length_string,
        "CN": f'{cnv["cn"]:g}',
        "LOG2": cnv["log2"],
        "DEPTH": cnv["depth"],
        "TYPE": cnv["type"],
        "TIER": tier,
        "REASON": reason,
        "GENES": ",".join(x["gene"] for x in reported_genes) or ".",
        "CDS_FRACTIONS": ",".join(f'{x["fraction"]:.4f}' for x in reported_genes) or ".",
        "DOSAGE_SCORES": ",".join("." if x["score"] is None else f'{x["score"]:.6g}' for x in reported_genes) or ".",
        "CLINVAR_ID": best_clinvar["id"] if best_clinvar else ".",
        "CLINVAR_PATIENT_FRACTION": f'{best_clinvar["patient_fraction"]:.4f}' if best_clinvar else ".",
        "CLINVAR_REFERENCE_FRACTION": f'{best_clinvar["reference_fraction"]:.4f}' if best_clinvar else ".",
        "CLINVAR_PHENOTYPE": best_clinvar["phenotype"] if best_clinvar else ".",
        "CLINVAR_REVIEW": best_clinvar["review"] if best_clinvar else ".",
        "DGV_ID": best_dgv["id"] if best_dgv else ".",
        "DGV_AF": "." if not best_dgv or best_dgv["af"] is None else f'{best_dgv["af"]:.6g}',
        "DGV_PATIENT_FRACTION": f'{best_dgv["patient_fraction"]:.4f}' if best_dgv else ".",
    })

    counts[tier] += 1

rows.sort(key=lambda x: (x["TIER"], chrom_key(x["CHROM"]), x["START"], x["END"]))

columns = [
    "CHROM", "START", "END", "COORDS", "LENGTH", "CN", "LOG2", "DEPTH", "TYPE",
    "TIER", "REASON", "GENES", "CDS_FRACTIONS", "DOSAGE_SCORES",
    "CLINVAR_ID", "CLINVAR_PATIENT_FRACTION", "CLINVAR_REFERENCE_FRACTION",
    "CLINVAR_PHENOTYPE", "CLINVAR_REVIEW",
    "DGV_ID", "DGV_AF", "DGV_PATIENT_FRACTION",
]

with open(args.output, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

print(f"CNVs passing basic filters: {len(cnvs):,}", file=sys.stderr)
print(f"Excluded by common DGV containment: {dgv_excluded:,}", file=sys.stderr)
print(f"Tier 1: {counts[1]:,}", file=sys.stderr)
print(f"Tier 2: {counts[2]:,}", file=sys.stderr)
print(f"Total prioritized CNVs: {len(rows):,}", file=sys.stderr)
