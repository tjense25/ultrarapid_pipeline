#!/usr/bin/env python3

import argparse
import base64
import csv
import html
import json
import re
from pathlib import Path

TABLE_RE = re.compile(r"^(\s*)const tableJson = (.+)$", re.MULTILINE)
SESSION_RE = re.compile(r"^(\s*)const sessionDictionary = (.+)$", re.MULTILINE)


def chrom_key(chrom):
    chrom = str(chrom).removeprefix("chr").upper()

    if chrom.isdigit():
        return int(chrom)

    return {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(chrom, 99)


def number(value, default=10**18):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def coordinates(values):
    chrom = values.get("CHROM", "")
    pos = values.get("POS", values.get("START", ""))

    if chrom:
        return chrom, number(pos)

    match = re.match(r"([^:]+):(\d+)", str(values.get("COORDS", "")))

    if match:
        return match.group(1), number(match.group(2))

    return "", 10**18


def read_igv_report(path):
    text = Path(path).read_text(encoding="utf-8")
    table_match = TABLE_RE.search(text)
    session_match = SESSION_RE.search(text)

    if not table_match or not session_match:
        raise SystemExit(f"Could not find embedded IGV data in {path}")

    table = json.loads(table_match.group(2))
    sessions = json.loads(session_match.group(2))

    return text, table, sessions


def merge_igv_reports(paths):
    template = None
    headers = None
    entries = []

    for path in paths:
        text, table, sessions = read_igv_report(path)

        if template is None:
            template = text
            headers = table["headers"]
        elif table["headers"] != headers:
            raise SystemExit(f"IGV report headers differ in {path}")

        for row in table["rows"]:
            old_id = str(row[0])

            if old_id not in sessions:
                raise SystemExit(f"Missing IGV session {old_id} in {path}")

            entries.append((list(row), sessions[old_id]))

    def sort_key(entry):
        values = dict(zip(headers, entry[0]))
        chrom, pos = coordinates(values)
        return number(values.get("TIER")), chrom_key(chrom), pos

    entries.sort(key=sort_key)

    rows = []
    sessions = {}

    for new_id, (row, session) in enumerate(entries):
        row[0] = new_id
        rows.append(row)
        sessions[str(new_id)] = session

    table = {
        "headers": headers,
        "rows": rows,
    }

    template = TABLE_RE.sub(
        lambda match: (
            f"{match.group(1)}const tableJson = "
            f"{json.dumps(table, separators=(',', ':'))}"
        ),
        template,
        count=1,
    )

    template = SESSION_RE.sub(
        lambda match: (
            f"{match.group(1)}const sessionDictionary = "
            f"{json.dumps(sessions, separators=(',', ':'))}"
        ),
        template,
        count=1,
    )

    return template.encode("utf-8")


def merge_tsvs(paths, output):
    fieldnames = None
    rows = []

    for path in paths:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")

            if fieldnames is None:
                fieldnames = reader.fieldnames
            elif reader.fieldnames != fieldnames:
                raise SystemExit(f"TSV headers differ in {path}")

            rows.extend(reader)

    if not fieldnames:
        raise SystemExit("No TSV headers found")

    def sort_key(row):
        chrom, pos = coordinates(row)
        return number(row.get("TIER")), chrom_key(chrom), pos

    rows.sort(key=sort_key)

    Path(output).parent.mkdir(parents=True, exist_ok=True)

    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


parser = argparse.ArgumentParser(
    description="Combine chunked IGV reports into one tabbed HTML report"
)
parser.add_argument("--sample", required=True)
parser.add_argument("--snv-html", nargs="+", required=True)
parser.add_argument("--sv-html", nargs="+", required=True)
parser.add_argument("--snv-tsv", nargs="+", required=True)
parser.add_argument("--sv-tsv", nargs="+", required=True)
parser.add_argument("--cnv-html", required=True)
parser.add_argument("--snv-tsv-output", required=True)
parser.add_argument("--sv-tsv-output", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

sample = html.escape(args.sample)

inputs = (
    args.snv_html
    + args.sv_html
    + args.snv_tsv
    + args.sv_tsv
    + [args.cnv_html]
)

for path in inputs:
    if not Path(path).is_file():
        raise SystemExit(f"Input not found: {path}")

merge_tsvs(args.snv_tsv, args.snv_tsv_output)
merge_tsvs(args.sv_tsv, args.sv_tsv_output)

reports = [
    ("snv", "SNV / Indel", merge_igv_reports(args.snv_html)),
    ("sv", "Structural Variants", merge_igv_reports(args.sv_html)),
    ("cnv", "Copy-Number Variants", Path(args.cnv_html).read_bytes()),
]

buttons = "\n".join(
    f'<button class="tab" data-report="{key}">{label}</button>'
    for key, label, _ in reports
)

buttons += '\n<button class="tab" data-report="methods">Filtering &amp; Tiers</button>'

panels = "\n".join(
    f'<section id="panel-{key}" class="panel">'
    f'<iframe id="frame-{key}" title="{label}"></iframe>'
    f"</section>"
    for key, label, _ in reports
)

payloads = "\n".join(
    f'<script type="application/octet-stream" id="payload-{key}">'
    f'{base64.b64encode(report).decode("ascii")}'
    f"</script>"
    for key, _, report in reports
)

methods = """
<section id="panel-methods" class="panel methods-panel">
<div class="methods-content">

<h2>Candidate filtering and tier definitions</h2>

<div class="method-card">
<h3>SNVs and small indels</h3>

<h4>General filtering</h4>
<ul>
    <li>Restricted to genes in the union of PanelApp Green and GenCC Strong disease-gene lists.</li>
    <li>VEP MODIFIER annotations are excluded.</li>
    <li>Variants classified as ClinVar Benign or Likely benign are excluded.</li>
    <li>SNV/indels with maximum gnomAD population allele frequency greater than 1% are excluded</li>
    <li>Indels overlapping GIAB homopolymer regions are excluded when QUAL is below 20.</li>
    <li>Other variants with QUAL below 10 are retained but flagged as LOW_QUAL.</li>
</ul>

<div class="tier tier1">
<span>Tier 1</span>
ClinVar Pathogenic or Likely pathogenic variants.
</div>

<div class="tier tier2">
<span>Tier 2</span>
High-impact variants in a selected disease gene that are absent from gnomAD
or have maximum population allele frequency below 0.001.
</div>
</div>

<div class="method-card">
<h3>Structural variants</h3>

<h4>General filtering</h4>
<ul>
    <li>Unresolved BND and TRA records are excluded.</li>
    <li>Population frequency for SVs evaluated using gnomAD-SV and the UWONT500 long-read reference.</li>
    <li>Gene overlap is restricted to supplied gene list (by default GenCC Definitive and Strong evidence genes).</li>
    <li>MANE coding overlap is determined using MANE Select coding-sequence intervals.</li>
</ul>

<div class="tier tier1">
<span>Tier 1</span>
Type-compatible overlap with a ClinVar Pathogenic or Likely pathogenic structural variant.
</div>

<div class="tier tier2">
<span>Tier 2</span>
Structural variants overlapping MANE coding sequence in a selected disease gene,
with gnomAD AF below 0.0001 or absent and no matching variant in UWONT500.
</div>
</div>

<div class="method-card">
<h3>Copy-number variants</h3>

<h4>General filtering</h4>
<ul>
    <li>Restricted to autosomal CNVkit segments at least 50 kb long.</li>
    <li>Calls must contain at least 10 coverage bins.</li>
    <li>Losses and gains are matched separately when comparing against reference CNVs.</li>
    <li>Common DGV Gold Standard CNVs with allele frequency above 1% are excluded.</li>
</ul>

<div class="tier tier1">
<span>Tier 1</span>
Type-matched overlap with a ClinVar Pathogenic or Likely pathogenic CNV,
provided the event is not common in DGV.
</div>

<div class="tier tier2">
<span>Tier 2</span>
A deletion covering at least 84% of the MANE Select CDS of a gene with
pHaplo ≥0.86, or a duplication covering at least 84% of the MANE Select CDS
of a gene with pTriplo ≥0.94. Defined from Collins et al 2022 10.1016/j.cell.2022.06.036
</div>
</div>

<div class="disclaimer">
Candidate tiers support rapid review and prioritization. They are not clinical
classifications and require manual review of read evidence, genomic context,
inheritance, phenotype relevance, and supporting databases.
</div>

</div>
</section>
"""

document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{sample} — Ultrarapid Variant Report</title>

<style>
:root {{
    font-family: Inter, Arial, sans-serif;
    color: #172033;
}}

* {{
    box-sizing: border-box;
}}

html, body {{
    width: 100%;
    height: 100%;
    margin: 0;
    overflow: hidden;
    background: white;
}}

.report-header {{
    height: 68px;
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding: 0 22px;
    color: white;
    background: linear-gradient(110deg, #18324a, #245c53);
}}

.report-header h1 {{
    margin: 0;
    font-size: 21px;
    font-weight: 650;
    letter-spacing: 0.01em;
}}

.sample-name {{
    padding: 8px 14px;
    border: 1px solid rgba(255,255,255,0.35);
    border-radius: 20px;
    background: rgba(255,255,255,0.12);
    font-size: 14px;
    font-weight: 600;
}}

.tabs {{
    height: 50px;
    display: flex;
    gap: 4px;
    padding: 5px 12px 0;
    background: #f1f5f9;
    border-bottom: 1px solid #cbd5e1;
}}

.tab {{
    padding: 0 20px;
    border: 1px solid transparent;
    border-radius: 7px 7px 0 0;
    background: transparent;
    color: #475569;
    font-size: 14px;
    font-weight: 600;
    cursor: pointer;
}}

.tab:hover {{
    color: #0f172a;
    background: #e2e8f0;
}}

.tab.active {{
    color: #0f172a;
    background: white;
    border-color: #cbd5e1;
    border-bottom-color: white;
}}

.panel {{
    display: none;
    position: absolute;
    inset: 118px 0 0;
}}

.panel.active {{
    display: block;
}}

iframe {{
    width: 100%;
    height: 100%;
    border: 0;
}}

.methods-panel {{
    overflow-y: auto;
    background: #f8fafc;
}}

.methods-content {{
    width: min(1050px, calc(100% - 40px));
    margin: 26px auto 50px;
}}

.methods-content h2 {{
    margin: 0 0 20px;
    font-size: 24px;
}}

.method-card {{
    margin-bottom: 18px;
    padding: 22px 25px;
    border: 1px solid #d8e0e8;
    border-radius: 10px;
    background: white;
    box-shadow: 0 2px 7px rgba(15,23,42,0.05);
}}

.method-card h3 {{
    margin: 0 0 16px;
    color: #18324a;
    font-size: 19px;
}}

.method-card h4 {{
    margin: 12px 0 6px;
    font-size: 14px;
}}

.method-card ul {{
    margin: 6px 0 18px;
    padding-left: 22px;
    line-height: 1.55;
}}

.tier {{
    margin-top: 10px;
    padding: 12px 15px;
    border-left: 5px solid;
    border-radius: 5px;
    line-height: 1.45;
}}

.tier span {{
    display: inline-block;
    margin-right: 8px;
    font-weight: 750;
}}

.tier1 {{
    border-color: #b91c1c;
    background: #fef2f2;
}}

.tier2 {{
    border-color: #d97706;
    background: #fffbeb;
}}

.disclaimer {{
    padding: 16px 20px;
    border: 1px solid #94a3b8;
    border-radius: 8px;
    color: #475569;
    background: #f1f5f9;
    font-size: 13px;
    line-height: 1.5;
}}
</style>
</head>

<body>

<header class="report-header">
    <h1>Ultrarapid Genome Variant Report</h1>
    <div class="sample-name">Sample: {sample}</div>
</header>

<nav class="tabs">
{buttons}
</nav>

{panels}
{methods}
{payloads}

<script>
function decodeReport(encoded) {{
    const chunkSize = 1024 * 1024;
    const chunks = [];

    for (let offset = 0; offset < encoded.length; offset += chunkSize) {{
        const binary = atob(encoded.slice(offset, offset + chunkSize));
        const bytes = new Uint8Array(binary.length);

        for (let i = 0; i < binary.length; i++) {{
            bytes[i] = binary.charCodeAt(i);
        }}

        chunks.push(bytes);
    }}

    return URL.createObjectURL(
        new Blob(chunks, {{type: "text/html;charset=utf-8"}})
    );
}}

function loadReport(name) {{
    if (name === "methods") return;

    const frame = document.getElementById(`frame-${{name}}`);

    if (frame.dataset.loaded) return;

    const payload = document.getElementById(`payload-${{name}}`);
    frame.src = decodeReport(payload.textContent.trim());
    frame.dataset.loaded = "true";

    payload.remove();
}}

function showReport(name) {{
    document.querySelectorAll(".tab").forEach(button => {{
        button.classList.toggle("active", button.dataset.report === name);
    }});

    document.querySelectorAll(".panel").forEach(panel => {{
        panel.classList.toggle("active", panel.id === `panel-${{name}}`);
    }});

    loadReport(name);
}}

document.querySelectorAll(".tab").forEach(button => {{
    button.addEventListener("click", () => showReport(button.dataset.report));
}});

showReport("snv");
</script>

</body>
</html>
"""

Path(args.output).parent.mkdir(parents=True, exist_ok=True)
Path(args.output).write_text(document, encoding="utf-8")

print(f"Wrote combined report: {args.output}")
print(f"Wrote combined SNV/indel table: {args.snv_tsv_output}")
print(f"Wrote combined SV table: {args.sv_tsv_output}")
