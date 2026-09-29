#!/usr/bin/env python3
import argparse
import base64
import html
from pathlib import Path

parser = argparse.ArgumentParser(description="Combine three IGV reports into one tabbed HTML report")
parser.add_argument("--sample", required=True)
parser.add_argument("--snv", required=True)
parser.add_argument("--sv", required=True)
parser.add_argument("--cnv", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

sample = html.escape(args.sample)

reports = [
    ("snv", "SNV / Indel", Path(args.snv)),
    ("sv", "Structural Variants", Path(args.sv)),
    ("cnv", "Copy-Number Variants", Path(args.cnv)),
]

for _, _, path in reports:
    if not path.is_file():
        raise SystemExit(f"Report not found: {path}")

buttons = "\n".join(
    f'<button class="tab" data-report="{key}">{label}</button>'
    for key, label, _ in reports
)

buttons += '\n<button class="tab" data-report="methods">Filtering &amp; Tiers</button>'

panels = "\n".join(
    f'<section id="panel-{key}" class="panel">'
    f'<iframe id="frame-{key}" title="{label}"></iframe>'
    f'</section>'
    for key, label, _ in reports
)

payloads = "\n".join(
    f'<script type="application/octet-stream" id="payload-{key}">'
    f'{base64.b64encode(path.read_bytes()).decode("ascii")}'
    f'</script>'
    for key, _, path in reports
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
of a gene with pTriplo ≥0.94. Defined from Collins et al 2022  10.1016/j.cell.2022.06.036
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

Path(args.output).write_text(document, encoding="utf-8")
print(f"Wrote combined report: {args.output}")
