#!/usr/bin/env python3

import argparse
import csv
import gzip
import re
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--genes", required=True)
parser.add_argument("--sample")
parser.add_argument("--common-af", type=float, default=0.01)
parser.add_argument("--ultrarare-af", type=float, default=0.001)
args = parser.parse_args()

columns = [
    "CHROM", "POS", "REF", "ALT", "RSID", "VARIANT_TYPE", "QUAL",
    "HOMOPOLYMER", "QUALITY_FLAG", "GT", "GQ", "DP", "VAF", "TIER",
    "REASON", "MAX_GNOMAD_AF", "GENE", "GENE_ID", "IMPACT",
    "CONSEQUENCE", "CLIN_SIG", "MIM_MORBID_DISEASE"
]

def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def numbers(value):
    output = []
    for token in re.split(r"[&,]", value or ""):
        try:
            output.append(float(token))
        except ValueError:
            pass
    return output

def max_gnomad_af(entry, gnomad_fields):
    values = [x for field in gnomad_fields for x in numbers(entry.get(field, ""))]
    return max(values) if values else None

def clinvar_terms(value):
    value = value.lower().replace(" ", "_").replace("-", "_")
    return {x for x in re.split(r"[&,/|]", value) if x}

def clinvar_benign(value):
    return bool(clinvar_terms(value) & {"benign", "likely_benign"})

def clinvar_plp(value):
    terms = clinvar_terms(value)
    excluded = {"benign", "likely_benign"}
    return bool(terms & {"pathogenic", "likely_pathogenic"}) and not bool(terms & excluded)

def clinvar_conflicting(value):
    terms = clinvar_terms(value)
    return bool(terms & {"uncertain_significance", "conflicting_classifications_of_pathogenicity", "conflicting_interpretations_of_pathogenicity"})

def transcript_rank(entry):
    if entry.get("MANE_SELECT") not in {"", ".", None}: return 0
    if entry.get("MANE_PLUS_CLINICAL") not in {"", ".", None}: return 1
    if entry.get("CANONICAL") == "YES": return 2
    return 3

def mim_morbid_diseases(entries):
    diseases = []
    for entry in entries:
        for annotation in entry.get("PHENOTYPES", "").split("&"):
            if not annotation:
                continue
            if annotation.startswith("MIM_morbid+"):
                parts = annotation.split("+")
                if len(parts) > 1:
                    diseases.append(parts[1])
            else:
                diseases.append(annotation)
    return ",".join(dict.fromkeys(diseases)) or "."

def chrom_key(chrom):
    value = chrom.removeprefix("chr")
    return int(value) if value.isdigit() else {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(value, 99)

with open(args.genes) as source:
    disease_genes = {line.strip().split()[0] for line in source if line.strip() and not line.startswith("#")}

disease_genes -= {"GENE", "Gene", "gene", "SYMBOL", "Symbol", "symbol"}

rows = []
counts = {1: 0, 2: 0}
input_records = 0
common_records = 0
benign_annotations = 0
modifier_annotations = 0
homopolymer_indels = 0
low_quality_records = 0
csq_fields = None
gnomad_fields = []
sample_index = None

with open_text(args.input) as source:
    for line in source:
        if line.startswith("##INFO=<ID=CSQ"):
            match = re.search(r'Format: ([^"]+)', line)
            if not match:
                raise SystemExit("Could not parse the VEP CSQ header")
            csq_fields = match.group(1).split("|")
            gnomad_fields = [x for x in csq_fields if x.startswith("gnomAD") and x.endswith("_AF")]
            print("Using gnomAD fields: " + ",".join(gnomad_fields), file=sys.stderr)
            continue

        if line.startswith("#CHROM"):
            header = line.rstrip().split("\t")
            if args.sample:
                if args.sample not in header:
                    raise SystemExit(f"Sample {args.sample} is absent from the VCF")
                sample_index = header.index(args.sample)
            else:
                sample_index = 9 if len(header) > 9 else None
            continue

        if line.startswith("#"):
            continue

        if csq_fields is None:
            raise SystemExit("No VEP CSQ header found")

        input_records += 1
        fields = line.rstrip().split("\t")
        if len(fields) < 8:
            continue

        chrom, pos, _, ref, alt, qual, record_filter, info = fields[:8]
        pos = int(pos)

        if record_filter != "PASS": continue
        info_fields = {}
        for item in info.split(";"):
            key, separator, value = item.partition("=")
            info_fields[key] = value if separator else "T"

        alleles = alt.split(",")
        is_indel = any(len(allele) != len(ref) for allele in alleles)
        variant_type = "INDEL" if is_indel else "SNV"
        homopolymer = "T" in re.split(r"[,&|]", info_fields.get("HOMOPOLYMER", ""))

        try:
            qual_number = float(qual)
        except ValueError:
            qual_number = None

        if is_indel and homopolymer and qual_number is not None and qual_number < 20:
            homopolymer_indels += 1
            continue

        quality_flag = "LOW_QUAL" if qual_number is not None and qual_number < 20 else "."
        if quality_flag == "LOW_QUAL":
            low_quality_records += 1

        csq = info_fields.get("CSQ")
        if not csq:
            continue

        entries = []
        for annotation in csq.split(","):
            values = annotation.split("|")
            values += [""] * (len(csq_fields) - len(values))
            entries.append(dict(zip(csq_fields, values)))

        all_af = [max_gnomad_af(entry, gnomad_fields) for entry in entries]
        all_af = [value for value in all_af if value is not None]
        record_max_af = max(all_af) if all_af else None

        if record_max_af is not None and record_max_af > args.common_af:
            common_records += 1
            continue

        sample = {}
        if sample_index is not None and len(fields) > sample_index:
            format_fields = fields[8].split(":")
            sample_values = fields[sample_index].split(":")
            sample = dict(zip(format_fields, sample_values))

        vaf = sample.get("AF", ".")
        best_by_gene = {}

        for entry in entries:
            gene = entry.get("SYMBOL", "")
            impact = entry.get("IMPACT", "").upper()
            clinical = entry.get("CLIN_SIG", "")

            if gene not in disease_genes:
                continue

            if impact == "MODIFIER":
                modifier_annotations += 1
                continue

            if clinvar_benign(clinical):
                benign_annotations += 1
                continue

            if clinvar_plp(clinical):
                tier = 1
                reason = "CLINVAR_PLP_CONFLICTING" if clinvar_conflicting(clinical) else "CLINVAR_PLP"
            elif impact == "HIGH" and (record_max_af is None or record_max_af < args.ultrarare_af):
                tier, reason = 2, "RARE_pLoF_DISEASE_GENE"
            else:
                continue

            candidate = (tier, transcript_rank(entry), entry, reason)
            if gene not in best_by_gene or candidate[:2] < best_by_gene[gene][:2]:
                best_by_gene[gene] = candidate

        for gene, (tier, _, entry, reason) in best_by_gene.items():
            gene_entries = [x for x in entries if x.get("SYMBOL") == gene]

            rows.append({
                "CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt, 
	 "RSID": ",".join(x for x in entry.get("Existing_variation", "").split("&") if x.startswith("rs")) or ".", 
	 "VARIANT_TYPE": variant_type,
                "QUAL": qual, "HOMOPOLYMER": "T" if homopolymer else ".", "QUALITY_FLAG": quality_flag, 
	 "GT": sample.get("GT", "."), "GQ": sample.get("GQ", "."), "DP": sample.get("DP", "."),
                "VAF": vaf, "TIER": tier, "REASON": reason, "MAX_GNOMAD_AF": "." if record_max_af is None else f"{record_max_af:.8g}",
                "GENE": gene, "GENE_ID": entry.get("Gene", ".") or ".", "IMPACT": entry.get("IMPACT", ".") or ".",
                "CONSEQUENCE": entry.get("Consequence", ".") or ".", "CLIN_SIG": entry.get("CLIN_SIG", ".") or ".",
                "MIM_MORBID_DISEASE": mim_morbid_diseases(gene_entries),
            })
            counts[tier] += 1

rows.sort(key=lambda x: (x["TIER"], chrom_key(x["CHROM"]), x["POS"]))

with open(args.output, "w", newline="") as destination:
    writer = csv.DictWriter(destination, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)

print(f"Input variants: {input_records:,}", file=sys.stderr)
print(f"Excluded common variants: {common_records:,}", file=sys.stderr)
print(f"Excluded ClinVar-benign annotations: {benign_annotations:,}", file=sys.stderr)
print(f"Excluded MODIFIER annotations: {modifier_annotations:,}", file=sys.stderr)
print(f"Excluded QUAL<20 homopolymer indels: {homopolymer_indels:,}", file=sys.stderr)
print(f"Retained LOW_QUAL variants: {low_quality_records:,}", file=sys.stderr)
print(f"Tier 1 candidates: {counts[1]:,}", file=sys.stderr)
print(f"Tier 2 candidates: {counts[2]:,}", file=sys.stderr)
