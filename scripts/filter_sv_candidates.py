#!/usr/bin/env python3
import argparse
import gzip
import os
import re
import subprocess
import tempfile
from urllib.parse import unquote_plus


p = argparse.ArgumentParser()
p.add_argument("vcf")
p.add_argument("genes")
p.add_argument("clinvar")
p.add_argument("output")
p.add_argument("--max-af", type=float, default=0.0001)
p.add_argument("--reciprocal-overlap", type=float, default=0.5)
a = p.parse_args()


def open_text(path):
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def parse_info(text):
    output = {}
    for item in text.split(";"):
        key, separator, value = item.partition("=")
        output[key] = value if separator else True
    return output


def numeric_values(value):
    output = []
    for token in re.split(r"[,|&]", str(value or "")):
        try:
            output.append(float(token))
        except ValueError:
            pass
    return output


def max_info(info, fields):
    values = []
    for field in fields:
        values.extend(numeric_values(info.get(field)))
    return max(values) if values else None


def unique_values(values):
    return list(dict.fromkeys(x for x in values if x and x != "."))


def csq_value(entries, field):
    values = unique_values(entry.get(field, "") for entry in entries)
    return ",".join(values) if values else "."


def summarize_genes(entries, field, limit=5):
    values = unique_values(entry.get(field, "") for entry in entries)
    if not values:
        return "."
    if len(values) > limit:
        return f"{len(values)} genes" if field == "SYMBOL" else f"{len(values)} gene IDs"
    return ",".join(values)

def mim_values(entries):
    values = []
    for entry in entries:
        for phenotype in entry.get("PHENOTYPES", "").split("&"):
            phenotype = phenotype.strip()
            if phenotype and phenotype != ".":
                values.append(unquote_plus(phenotype))
    values = unique_values(values)
    return ",".join(values) if values else "."

def clinical_type(text):
    text = text.lower().replace("_", " ")

    if "copy number loss" in text or "deletion" in text:
        return "DEL"
    if "copy number gain" in text or "duplication" in text:
        return "DUP"
    if "inversion" in text:
        return "INV"
    if "insertion" in text:
        return "INS"
    if "copy number" in text or "cnv" in text:
        return "CNV"

    return "UNK"


def compatible(query, clinical):
    return clinical == query or (
        clinical == "CNV" and query in {"DEL", "DUP", "CNV"}
    )


def clinvar_plp(text):
    text = text.lower().replace("_", " ")

    return (
        "pathogenic" in text
        and "benign" not in text
        and "uncertain" not in text
        and "conflict" not in text
    )


def chrom_key(chrom):
    value = chrom.removeprefix("chr")

    if value.isdigit():
        return int(value)

    return {
        "X": 23,
        "Y": 24,
        "M": 25,
        "MT": 25,
    }.get(value, 99)


def clinvar_description(fields, reciprocal_overlap):
    clinvar_id = fields[3] or "."
    clinical_significance = fields[9] or "."
    variant_type = fields[10] or "."
    phenotype = fields[18] or fields[12] or "."
    review = fields[21] or "."

    phenotype = phenotype.replace("|", "/")
    review = review.replace("|", "/")

    return (
        f"{clinvar_id}"
        f"|{variant_type}"
        f"|{clinical_significance}"
        f"|reciprocal_overlap={reciprocal_overlap:.4f}"
        f"|{phenotype}"
        f"|{review}"
    )


genes = set()

with open(a.genes) as source:
    for line in source:
        if line.strip() and not line.startswith("#"):
            genes.add(line.split()[0])

genes -= {
    "gene",
    "Gene",
    "GENE",
    "symbol",
    "Symbol",
    "SYMBOL",
}


records = []
csq_fields = None
query_bed = None

try:
    with tempfile.NamedTemporaryFile(
        "w",
        suffix=".bed",
        delete=False,
    ) as query:
        query_bed = query.name

        with open_text(a.vcf) as source:
            for line in source:
                if line.startswith("##INFO=<ID=CSQ"):
                    match = re.search(r'Format: ([^"]+)', line)
                    csq_fields = match.group(1).split("|") if match else None
                    continue

                if line.startswith("#"):
                    continue

                if not csq_fields:
                    raise SystemExit("Could not find VEP CSQ header")

                columns = line.rstrip().split("\t")

                if len(columns) < 8:
                    continue

                info = parse_info(columns[7])
                svtype = info.get("SVTYPE", "").upper()

                if (
                    svtype in {"BND", "TRA"}
                    or "[" in columns[4]
                    or "]" in columns[4]
                ):
                    continue

                csq = []

                for annotation in info.get("CSQ", "").split(","):
                    if not annotation:
                        continue

                    values = annotation.split("|")
                    values += [""] * (len(csq_fields) - len(values))
                    csq.append(dict(zip(csq_fields, values)))

                position = int(columns[1])
                end = int(float(info.get("END", position)))
                start0 = position - 1
                end0 = max(end, start0 + 1)
                index = len(records)

                records.append((columns, info, csq))

                query.write(
                    f"{columns[0]}\t{start0}\t{end0}\t"
                    f"{index}\t{svtype}\n"
                )

    hits = {}

    command = [
        "bedtools",
        "intersect",
        "-wa",
        "-wb",
        "-a",
        query_bed,
        "-b",
        a.clinvar,
    ]

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        text=True,
    )

    for line in process.stdout:
        fields = line.rstrip().split("\t")

        if len(fields) < 27:
            continue

        query_chrom = fields[0]
        query_start = int(fields[1])
        query_end = int(fields[2])
        query_index = int(fields[3])
        query_type = fields[4]

        clinvar = fields[5:]
        clinvar_chrom = clinvar[0]
        clinvar_start = int(clinvar[1])
        clinvar_end = int(clinvar[2])
        clinical_significance = clinvar[9]
        clinvar_type = clinical_type(clinvar[10])

        if query_chrom != clinvar_chrom:
            continue

        if not clinvar_plp(clinical_significance):
            continue

        if not compatible(query_type, clinvar_type):
            continue

        overlap = max(
            0,
            min(query_end, clinvar_end)
            - max(query_start, clinvar_start),
        )

        patient_fraction = overlap / max(1, query_end - query_start)
        clinvar_fraction = overlap / max(
            1,
            clinvar_end - clinvar_start,
        )
        reciprocal_overlap = min(
            patient_fraction,
            clinvar_fraction,
        )

        if reciprocal_overlap < a.reciprocal_overlap:
            continue

        hit = {
            "score": reciprocal_overlap,
            "description": clinvar_description(
                clinvar,
                reciprocal_overlap,
            ),
        }

        if (
            query_index not in hits
            or reciprocal_overlap > hits[query_index]["score"]
        ):
            hits[query_index] = hit

    if process.wait() != 0:
        raise SystemExit("bedtools intersect failed")

finally:
    if query_bed and os.path.exists(query_bed):
        os.unlink(query_bed)


header = [
    "CHROM",
    "POS",
    "END",
    "COORDS",
    "ID",
    "SVTYPE",
    "SVLEN",
    "QUAL",
    "FILTER",
    "GT",
    "gnomad_AF",
    "UWONT500_COUNT",
    "TIER",
    "GENE",
    "GENE_ID",
    "VEP_IMPACT",
    "VEP_CONSEQUENCE",
    "CLINVAR_SV",
    "MIM_MORBID"
]


rows = []
counts = {
    1: 0,
    2: 0,
}


for index, (columns, info, csq) in enumerate(records):
    matched = [
        entry
        for entry in csq
        if entry.get("SYMBOL") in genes
    ]

    mane_cds = "T" in re.split(
        r"[,|&]",
        str(info.get("MANE_CDS", "")),
    )

    max_af = max_info(
        info,
        ("Max_AF", "Max_PopMax_AF"),
    )

    uwont_count = max_info(
        info,
        ("UWONT500_Count",),
    )

    if uwont_count is None:
        uwont_count = 0

    if index in hits:
        tier = 1
        selected = csq
    elif (
        matched
        and mane_cds
        and (max_af is None or max_af < a.max_af)
        and uwont_count < 1
    ):
        tier = 2
        selected = matched
    else:
        continue

    original_position = int(columns[1])
    original_end = int(float(info.get("END", original_position)))

    svlen_values = numeric_values(info.get("SVLEN"))

    if svlen_values:
        svlen = int(svlen_values[0])
    else:
        svlen = original_end - original_position + 1

    coords = (
        f"{columns[0]}:"
        f"{original_position}-"
        f"{original_end}"
    )

    report_start = max(0, original_position - 1000)

    if abs(svlen) < 50000:
        report_end = original_end + 1000
    else:
        report_end = original_position + 1000

    gt = dict(zip(columns[8].split(":"), columns[9].split(":"))).get("GT", ".") if len(columns) > 9 else "."

    rows.append([
        columns[0],
        str(report_start),
        str(report_end),
        coords,
        columns[2],
        info.get("SVTYPE", "."),
        info.get("SVLEN", "."),
        columns[5],
        columns[6],
        gt,
        "." if max_af is None else f"{max_af:g}",
        f"{uwont_count:g}",
        str(tier),
        summarize_genes(selected, "SYMBOL"),
        summarize_genes(selected, "Gene"),
        csq_value(selected, "IMPACT"),
        csq_value(selected, "Consequence"),
        hits[index]["description"] if index in hits else ".",
        mim_values(csq)
    ])

    counts[tier] += 1


rows.sort(
    key=lambda row: (
        int(row[12]),
        chrom_key(row[0]),
        int(row[1]),
    )
)


with open(a.output, "w") as output:
    output.write("\t".join(header) + "\n")

    for row in rows:
        output.write("\t".join(row) + "\n")


print(f"Tier 1: {counts[1]:,}")
print(f"Tier 2: {counts[2]:,}")
