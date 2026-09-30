import glob
from pathlib import Path

SAMPLE = config["sample"]
BAM_FILES=[]
if "bam_dir" in config: 
    BAM_FILES = [str(path) for path in Path(config["bam_dir"]).rglob("*.bam")]
    print("Merging %d bam files in %s" % (len(BAM_FILES), config["bam_dir"]))
INPUT_BAM=config["bam"] if "bam" in config else SAMPLE+"/merged_bam/"+SAMPLE+".merged.bam"
REFERENCE=config["fasta"]
PORE_VERSION=config["pore_version"].lower()
BASECALL_VERSION=config["basecall_version"].lower()
CLAIR3_MODELS = {
    ("r9", "hac"): "r941_prom_hac_g360+g422",
    ("r9", "sup"): "r941_prom_sup_g5014",
    ("r10", "hac"): "r1041_e82_400bps_hac_v500",
    ("r10", "sup"): "r1041_e82_400bps_sup_v500",
}
ANNOT_DIR=config["annotation_dir"]
GENE_LIST=config["gene_list"] if ("gene_list" in config and Path(config["gene_list"]).is_file()) else ANNOT_DIR+"/GenCC_DefinitiveStrong_genelist.txt"
MODEL=CLAIR3_MODELS[(PORE_VERSION,BASECALL_VERSION)]
gpus=int(config["gpus"])
print("Using model %s for Clair3 snv/indel calling" % MODEL)
print("Using %d GPUs for clair3 variant calling" % gpus)
print("Reading variant annotations from %s" % ANNOT_DIR)
print("Using %s gene list for variant prioritization" % GENE_LIST)

rule all:
    input:
        expand("{sample}/candidates/{sample}.cnv.candidate_variants.igv_report.no_snvs.html",sample=SAMPLE),
        expand("{sample}/candidates/{sample}.all_candidate_variants.igv_report.html",sample=SAMPLE)

rule merge_bams:
    input:
        BAM_FILES
    output:
        bam="{sample}/merged_bam/{sample}.merged.bam",
        bai="{sample}/merged_bam/{sample}.merged.bam.bai"
    threads: 128
    conda: "envs/ultrarapid.yaml"
    benchmark: "{sample}/benchmarks/merge.tsv"
    shell:
        "samtools merge -@ {threads} -o {output.bam} {input} && samtools index -@ {threads} {output.bam}"

rule clair3:
    input:
        bam=INPUT_BAM,
        ref=REFERENCE
    params:
        model=ANNOT_DIR+"/clair3_models/" + MODEL,
        gpu_args="--use_gpu --device=cuda:"+','.join(map(str,range(gpus))) if gpus > 0 else "",
        outdir="{sample}/clair3_out"
    output:
        "{sample}/clair3_out/{sample}.merge_output.PASS.vcf.gz"
    threads: 128
    benchmark: "{sample}/benchmarks/clair3.tsv"
    conda: "envs/ultrarapid.yaml"
    shell: """
        run_clair3.sh --bam_fn={input.bam} --ref_fn={input.ref} \
          --threads={threads} --platform=ont --model_path={params.model} \
          --output={params.outdir} --use_longphase_for_intermediate_phasing --min_coverage 4 --qual 10 \
          {params.gpu_args}
        bcftools view --threads 32 -f PASS -Oz -o {output} {params.outdir}/merge_output.vcf.gz
        bcftools index --threads 32 -t {output}
    """

rule sniffles:
    input:
        INPUT_BAM
    params:
        raw="{sample}/sniffles/{sample}.sniffles.vcf"
    output:
        "{sample}/sniffles/{sample}.sniffles.vcf.gz"
    threads: 16
    benchmark: "{sample}/benchmarks/sniffles.tsv"
    conda: "envs/ultrarapid.yaml"
    shell: """ 
        sniffles --input {input} --vcf {params.raw} --threads {threads} --minsupport 2 --minsvlen 30 --output-rnames --allow-overwrite
        bcftools sort -Oz -o {output} {params.raw}
        bcftools index --threads {threads} --tbi {output}
        rm -f {params.raw}
    """

rule cnvkit:
    input:
        INPUT_BAM
    params:
        flat_ref=ANNOT_DIR+"/GRCh38.cnvkit.flat_reference.cnn",
        outdir="{sample}/cnvkit",
        tmp_ref="{sample}/cnvkit/tmp.cnvkit.ref_regions.bed",
        cnn="{sample}/cnvkit/{sample}.cnn",
        cnr="{sample}/cnvkit/{sample}.cnr",
        cns="{sample}/cnvkit/{sample}.cns",
        call="{sample}/cnvkit/{sample}.cns.called"
    output:
        cnv="{sample}/cnvkit/{sample}.cnvkit.germline_autosomal_CNVs.called.bed",
        cov_bars="{sample}/cnvkit/{sample}.5kb_coverage_bars.bedgraph.gz"
    threads: 32
    benchmark: "{sample}/benchmarks/cnvkit.tsv"
    conda: "envs/ultrarapid.yaml"
    shell: """
        #get per base and per 5kb region depth quickly with mosdepth
        mosdepth -t {threads} --fast-mode --by 5000 {params.outdir}/mosdepth {input}
        tabix -p bed {params.outdir}/mosdepth.per-base.bed.gz

        zcat {params.outdir}/mosdepth.regions.bed.gz | awk 'BEGIN{{OFS="\t"}}{{print $1,$2,$3,$4}}' | bgzip -@ {threads} -c > {output.cov_bars}
        tabix -f -p bed {output.cov_bars}
 
        sed '1d' {params.flat_ref} > {params.tmp_ref}
        cnvkit.py coverage {params.outdir}/mosdepth.per-base.bed.gz {params.tmp_ref} --output {params.cnn}
        cnvkit.py fix {params.cnn} -r {params.flat_ref} --output {params.cnr}
        cnvkit.py segment {params.cnr} --method hmm-germline --output {params.cns}
        cnvkit.py call {params.cns} --method clonal --ploidy 2 --output {params.call}
        awk 'BEGIN{{FS=OFS="\t"}} NR==1 || (($3-$2)>=50000 && $6!=2 && $8>=10 && (($6<2 && $5<=-0.6) || ($6>2 && $5>=0.4)))' {params.call} | \
	 sed 's/^chrom/#chrom/' > {output.cnv} #make it valid bed by commenting out header
        rm {params.tmp_ref}
    """

rule longphase:
    input:
        snv="{sample}/clair3_out/{sample}.merge_output.PASS.vcf.gz",
        sv="{sample}/sniffles/{sample}.sniffles.vcf.gz",
        bam=INPUT_BAM,
        ref=REFERENCE
    params:
        prefix="{sample}/longphase/{sample}.phased"
    output:
        snv="{sample}/longphase/{sample}.phased.vcf.gz",
        sv="{sample}/longphase/{sample}.phased_SV.vcf.gz"
    threads: 64
    benchmark: "{sample}/benchmarks/longphase.tsv"
    conda: "envs/ultrarapid.yaml"
    shell: """ 
        longphase phase -s {input.snv} --sv-file {input.sv} -b {input.bam} -r {input.ref} -t {threads} -o {params.prefix} --ont
        bcftools sort -m 4G -Oz -o {output.snv} {params.prefix}.vcf
        bcftools index --threads 32 --tbi {output.snv}
        bcftools sort -m 4G -Oz -o {output.sv} {params.prefix}_SV.vcf
        bcftools index --threads 32 --tbi {output.sv}
        rm -f {params.prefix}*.vcf 
    """

rule haplotag:
    input:
        snv="{sample}/longphase/{sample}.phased.vcf.gz",
        sv="{sample}/longphase/{sample}.phased_SV.vcf.gz",
        bam=INPUT_BAM,
        ref=REFERENCE
    params:
        prefix="{sample}/longphase/{sample}.haplotagged"
    output:
        "{sample}/longphase/{sample}.haplotagged.bam"
    threads: 48
    benchmark: "{sample}/benchmarks/haplotag.tsv"
    conda: "envs/ultrarapid.yaml"
    shell: 
        "longphase haplotag -s {input.snv}  --sv-file {input.sv} -b {input.bam}  -r {input.ref}  -t {threads} -o {params.prefix} --tagSupplementary && samtools index -@ {threads} {output}"


rule vep_snv:
    input:
        snv="{sample}/longphase/{sample}.phased.vcf.gz"     
    params:
        cache=ANNOT_DIR+"/VEP",
        homopolymers=ANNOT_DIR+"/GRCh38_all_homopolymers_ge4_tagged.bed.gz",
        tmp_vcf="{sample}/vep/tmp.{sample}.vep.vcf"
    output:
        "{sample}/vep/{sample}.phased.vep.homo_anno.vcf.gz"
    threads: 96
    benchmark: "{sample}/benchmarks/vep.tsv"
    conda: "envs/vep.yaml"
    shell: """ 
        vep --input_file {input} --output_file {params.tmp_vcf} --check_existing \
          --format vcf --vcf --offline --compress_output bgzip --species homo_sapiens \
          --assembly GRCh38 --cache --cache_version 116 --dir_cache {params.cache} \
          --dir_plugins {params.cache}/Plugins --fork {threads} --buffer_size 20000 --symbol \
          --mane --canonical --biotype --variant_class --pick_allele_gene \
          --pick_order mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,rank,length \
          --check_existing --clin_sig_allele 1 --af_gnomade --af_gnomadg --max_af \
          --plugin "Phenotypes,file={params.cache}/phenotypes_116_GRCh38.gff.gz,include_types=Gene,include_sources=MIM_morbid,phenotype_feature=1,id_match=1,cols=source&phenotype&id" \
          --no_stats --force_overwrite
        bcftools annotate --threads 32 -a  {params.homopolymers} -h <(printf '##INFO=<ID=HOMOPOLYMER,Number=1,Type=String,Description="Variant overlaps GIAB GRCh38 homopolymer region">\n') \
            -c CHROM,FROM,TO,INFO/HOMOPOLYMER -Oz -o {output} {params.tmp_vcf}
        bcftools index --threads {threads} --tbi {output}
        rm -f {params.tmp_vcf}
    """

rule sv_annotate:
    input:
        "{sample}/longphase/{sample}.phased_SV.vcf.gz"
    params:
        ont_ref=ANNOT_DIR + "/UWONT_500_ONT.svafotate.bed.gz",
        core_ref=ANNOT_DIR + "/SVAFotate_core_SV_popAFs.GRCh38.bed.gz",
        vep_cache=ANNOT_DIR + "/VEP",
        MANE_CDS=ANNOT_DIR + "/MANE_CDS.bed.gz",
        tmp_cds = "{sample}/sv_anno/tmp.{sample}.cds_overlap.vcf",
        tmp_ont = "{sample}/sv_anno/tmp.{sample}.cds_overlap.ontanno.vcf",
        tmp_afs = "{sample}/sv_anno/tmp.{sample}.cds_overlap.svafotate.vcf"
    output:
        "{sample}/sv_anno/{sample}.phased_SV.svafotate.vep.vcf.gz"
    threads: 16
    benchmark: "{sample}/benchmarks/svafotate.tsv"
    conda: "envs/vep.yaml"
    shell: """ 
        #annotate with ONT SVs from UWONT500
        bcftools annotate --threads {threads} -a {params.MANE_CDS} -h <(echo '##INFO=<ID=MANE_CDS,Number=1,Type=String,Description="Overlaps MANE CDS">') \
          -c CHROM,FROM,TO,INFO/MANE_CDS -Ov -o {params.tmp_cds} {input}
        svafotate annotate -v {params.tmp_cds} -b {params.ont_ref} -o {params.tmp_ont} -f 0.5 --cpu {threads}
        sed -i -e 's/Max_AF/ONT_AF/' -e 's/Max_Het/ONT_Het/' -e 's/Max_HomAlt/ONT_HomAlt/' {params.tmp_ont}
        #annotate with short-read SVs from gnomad/1000G/CCDG
        svafotate annotate -v {params.tmp_ont} -b {params.core_ref} -o {params.tmp_afs} -f 0.5 --cpu {threads} 
        
        #annotate VEP consequence and MIM morbid
        vep --input_file {params.tmp_afs} --output_file {output} \
            --vcf --compress_output bgzip --force_overwrite --offline --cache --dir_cache {params.vep_cache} \
            --dir_plugins {params.vep_cache}/Plugins --assembly GRCh38 --fork {threads} --buffer_size 10000 \
            --symbol --canonical --mane --biotype --variant_class --pick_allele_gene --pick_order mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,rank,length \
            --plugin Phenotypes,file={params.vep_cache}/phenotypes_116_GRCh38.gff.gz,include_types=Gene,include_sources=MIM_morbid
        bcftools index --threads {threads} --tbi {output}
        rm -f {params.tmp_cds} {params.tmp_ont} {params.tmp_afs}
    """

rule snvindel_report:
    input:
        vcf="{sample}/vep/{sample}.phased.vep.homo_anno.vcf.gz",
        bam="{sample}/longphase/{sample}.haplotagged.bam"
    params:
        gene_list=GENE_LIST,
        igv_json="{sample}/candidates/tmp.{sample}.snvindel.igv_track.json"
    output:
        tsv="{sample}/candidates/{sample}.snv_indel.candidate_variants.tsv",
        html="{sample}/candidates/{sample}.snv_indel.candidate_variants.igv_report.html"
    threads: 1
    benchmark: "{sample}/benchmarks/snvindel_report.tsv"
    conda: "envs/igvreport.yaml"
    shell: """ 
        python3 scripts/filter_snv_candidates.py \
                --input {input.vcf} \
                --output {output.tsv} \
                --genes {params.gene_list}
        echo '[
            {{"name":"Phased VEP variants","type":"variant","format":"vcf","url":"{input.vcf}","displayMode":"EXPANDED"}}, 
            {{"name":"Haplotagged reads","type":"alignment","format":"bam","url":"{input.bam}","groupBy":"tag:HP","colorBy":"strand","showCoverage":true,"showAlignments":true,"height":500}}
         ]' > {params.igv_json}
        create_report {output.tsv}  --genome hg38 --sequence 1 --begin 2 --end 2 --flanking 1000 \
             --track-config {params.igv_json} --title "Ultrarapid candidate variants" --output {output.html}
        rm -f {params.igv_json}
    """

rule sv_report:
    input:
        vcf="{sample}/sv_anno/{sample}.phased_SV.svafotate.vep.vcf.gz",
        cnv="{sample}/cnvkit/{sample}.cnvkit.germline_autosomal_CNVs.called.bed",
        bam="{sample}/longphase/{sample}.haplotagged.bam"
    params:
        gene_list=GENE_LIST,
        clinvar_svs=ANNOT_DIR + "/nstd102_clinvar_pathogenic.GRCh38.bed.gz",
        clinvar_sv_igv=ANNOT_DIR + "/ClinVar_pathogenic_CNVs.SVs_lt50kb.igv.bed.gz",
        clinvar_cnv_igv=ANNOT_DIR + "/ClinVar_pathogenic_CNVs.largeCNVs.igv.bed.gz",
        igv_json="{sample}/candidates/tmp.{sample}.sv.igv_track.json"
    output:
        tsv="{sample}/candidates/{sample}.sv.candidate_variants.tsv",
        html="{sample}/candidates/{sample}.sv.candidate_variants.igv_report.html"
    threads: 1
    benchmark: "{sample}/benchmarks/sv_report.tsv"
    conda: "envs/igvreport.yaml"
    shell: """ 
        python3 scripts/filter_sv_candidates.py {input.vcf} {params.gene_list} {params.clinvar_svs} {output.tsv}
        echo '[
          {{"name":"Phased patient SVs","type":"variant","format":"vcf","url":"{input.vcf}","displayMode":"EXPANDED"}}, 
          {{"name":"CNVkit patient CNVs","type":"annotation","format":"bed","url": "{input.cnv}","displayMode":"SQUISHED"}},
          {{"name":"ClinVar pathogenic SVs","type":"annotation","format":"bed","url":"{params.clinvar_sv_igv}","displayMode":"SQUISHED","height":50}},
          {{"name":"ClinVar pathogenic CNVs","type":"annotation","format":"bed","url":"{params.clinvar_cnv_igv}","displayMode":"SQUISHED","height":100}},
          {{"name":"Haplotagged reads","type":"alignment","format":"bam","url":"{input.bam}","groupBy":"tag:HP","colorBy":"strand","showCoverage":true,"showAlignments":true,"height":600}}
         ]' > {params.igv_json}
        create_report {output.tsv}  --genome hg38 --sequence 1 --begin 2 --end 3 --flanking 5000 \
            --track-config {params.igv_json} --title "Ultrarapid candidate SVs" --output {output.html} \
            --info-columns COORDS ID SVTYPE SVLEN QUAL GT gnomad_AF UWONT500_COUNT TIER GENE GENE_ID VEP_IMPACT VEP_CONSEQUENCE CLINVAR_SV MIM_MORBID
        rm -f {params.igv_json}
    """

rule cnv_report_no_snvs:
    input:
        cnv="{sample}/cnvkit/{sample}.cnvkit.germline_autosomal_CNVs.called.bed",
        cov_bars="{sample}/cnvkit/{sample}.5kb_coverage_bars.bedgraph.gz",
        sniffles="{sample}/sniffles/{sample}.sniffles.vcf.gz"
    params:
        gene_list=GENE_LIST,
        clinvar_cnvs=ANNOT_DIR + "/nstd102_clinvar_pathogenic.GRCh38.bed.gz",
        dgv_cnvs=ANNOT_DIR + "/DGV_GoldStandard_hg38.bed.gz",
        collins_tsv=ANNOT_DIR + "/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz",
        mane_gtf=ANNOT_DIR + "/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz",
        clinvar_cnv_igv=ANNOT_DIR + "/ClinVar_pathogenic_CNVs.largeCNVs.igv.bed.gz",
        dgv_igv=ANNOT_DIR + "/DGV_GoldStandard_CNVs.igv.bed.gz",
        collins_igv=ANNOT_DIR + "/Collins_dosage_sensitivity.GRCh38.bed.gz",
        igv_json="{sample}/candidates/tmp.{sample}.cnvs.igv_track.json"
    output:
        tsv="{sample}/candidates/{sample}.cnv.candidate_variants.tsv",
        html="{sample}/candidates/{sample}.cnv.candidate_variants.igv_report.no_snvs.html"
    threads: 1
    benchmark: "{sample}/benchmarks/cnvreport_w_snvs.txt"
    conda: "envs/igvreport.yaml"
    shell: """
        python scripts/filter_cnv_candidates.py --cnv {input.cnv} \
	 --clinvar {params.clinvar_cnvs} --dgv {params.dgv_cnvs} \
	 --collins {params.collins_tsv} --mane {params.mane_gtf} \
	 --genes {params.gene_list} --output {output.tsv}

        MEDIAN=$(zcat {input.cov_bars} | awk '$1!="chrX" && $1!="chrY" && $1!="chrM" && $1!="chrMT"{{print $4}}' | sort -n | awk '{{x[NR]=$1}} END{{print NR%2 ? x[(NR+1)/2] : (x[NR/2]+x[NR/2+1])/2}}') 
        YMAX=$(awk -v m="$MEDIAN" 'BEGIN{{print 10*int((2.2*m+9.999999)/10)}}')
        echo "Median=$MEDIAN  IGV_max=$YMAX"
        echo '[
          {{"name":"CNVkit patient CNVs","type":"annotation","format":"bed","url": "{input.cnv}","displayMode":"EXPANDED","color": "rgb(34,139,94)"}},
          {{"name":"Mean read depth (5kb bins)","type":"wig","format":"bedgraph","url": "{input.cov_bars}","graphType":"bar","autoscale":false,"min":0,"max":'$YMAX',"color":"rgb(110,110,110)","height": 250}},
          {{"name":"Sniffles structural variants","type":"variant","format":"vcf","url":"{input.sniffles}","displayMode":"EXPANDED"}},
          {{"name":"ClinVar pathogenic CNVs","type":"annotation","format":"bed","url":"{params.clinvar_cnv_igv}","displayMode":"SQUISHED","height":200}},
          {{"name": "DGV Gold Standard CNVs","type":"annotation","format":"bed","url":"{params.dgv_igv}","displayMode":"SQUISHED","height":100}},
          {{"name": "Collins et al. dosage-sensitive genes","type":"annotation","format": "bed","url":"{params.collins_igv}","displayMode":"EXPANDED"}}
        ]' > {params.igv_json}
        create_report {output.tsv} --genome hg38 --sequence 1 --begin 2 --end 3 --flanking 1000000 \
	 --track-config {params.igv_json} --title "Ultrarapid candidate CNVs" --output {output.html} \
	 --info-columns COORDS LENGTH CN LOG2 DEPTH TYPE TIER GENES CDS_FRACTIONS DOSAGE_SCORES CLINVAR_ID CLINVAR_PATIENT_FRACTION CLINVAR_REFERENCE_FRACTION CLINVAR_PHENOTYPE CLINVAR_REVIEW DGV_ID DGV_AF DGV_PATIENT_FRACTION
        rm -f {params.igv_json}
    """

rule cnv_report:
    input:
        cnv="{sample}/cnvkit/{sample}.cnvkit.germline_autosomal_CNVs.called.bed",
        cand_tsv="{sample}/candidates/{sample}.cnv.candidate_variants.tsv",
        cov_bars="{sample}/cnvkit/{sample}.5kb_coverage_bars.bedgraph.gz",
        snvs="{sample}/clair3_out/{sample}.merge_output.PASS.vcf.gz",
        sniffles="{sample}/sniffles/{sample}.sniffles.vcf.gz"
    params:
        clinvar_cnv_igv=ANNOT_DIR + "/ClinVar_pathogenic_CNVs.largeCNVs.igv.bed.gz",
        dgv_igv=ANNOT_DIR + "/DGV_GoldStandard_CNVs.igv.bed.gz",
        collins_igv=ANNOT_DIR + "/Collins_dosage_sensitivity.GRCh38.bed.gz",
        igv_json="{sample}/candidates/tmp.{sample}.cnvs.igv_track.json"
    output:
        html="{sample}/candidates/{sample}.cnv.candidate_variants.igv_report.html"
    threads: 1
    benchmark: "{sample}/benchmarks/cnvreport.txt"
    conda: "envs/igvreport.yaml"
    shell: """
        MEDIAN=$(zcat {input.cov_bars} | awk '$1!="chrX" && $1!="chrY" && $1!="chrM" && $1!="chrMT"{{print $4}}' | sort -n | awk '{{x[NR]=$1}} END{{print NR%2 ? x[(NR+1)/2] : (x[NR/2]+x[NR/2+1])/2}}') 
        YMAX=$(awk -v m="$MEDIAN" 'BEGIN{{print 10*int((2.2*m+9.999999)/10)}}')
        echo "Median=$MEDIAN  IGV_max=$YMAX"
        echo '[
          {{"name":"CNVkit patient CNVs","type":"annotation","format":"bed","url": "{input.cnv}","displayMode":"EXPANDED","color":"rgb(34,139,94)"}},
          {{"name":"Mean read depth (5kb bins)","type":"wig","format":"bedgraph","url": "{input.cov_bars}","graphType":"bar","autoscale":false,"min":0,"max":'$YMAX',"color":"rgb(110,110,110)","height": 250}},
          {{"name":"Sniffles structural variants","type":"variant","format":"vcf","url":"{input.sniffles}","displayMode":"EXPANDED"}},
          {{"name":"Clair3 SNV/indels","type":"variant","format":"vcf","url":"{input.snvs}","displayMode":"EXPANDED"}}, 
          {{"name":"ClinVar pathogenic CNVs","type":"annotation","format":"bed","url":"{params.clinvar_cnv_igv}","displayMode":"SQUISHED","height":200}},
          {{"name": "DGV Gold Standard CNVs","type":"annotation","format":"bed","url":"{params.dgv_igv}","displayMode":"SQUISHED","height":100}},
          {{"name": "Collins et al. dosage-sensitive genes","type":"annotation","format": "bed","url":"{params.collins_igv}","displayMode":"EXPANDED"}}
        ]' > {params.igv_json}
        create_report {input.cand_tsv} --genome hg38 --sequence 1 --begin 2 --end 3 --flanking 1000000 \
	 --track-config {params.igv_json} --title "Ultrarapid candidate CNVs" --output {output.html} \
	 --info-columns COORDS LENGTH CN LOG2 DEPTH TYPE TIER GENES CDS_FRACTIONS DOSAGE_SCORES CLINVAR_ID CLINVAR_PATIENT_FRACTION CLINVAR_REFERENCE_FRACTION CLINVAR_PHENOTYPE CLINVAR_REVIEW DGV_ID DGV_AF DGV_PATIENT_FRACTION
        rm -f {params.igv_json}
    """

rule combine_igv_reports:
    input:
        snv="{sample}/candidates/{sample}.snv_indel.candidate_variants.igv_report.html",
        sv="{sample}/candidates/{sample}.sv.candidate_variants.igv_report.html",
        cnv="{sample}/candidates/{sample}.cnv.candidate_variants.igv_report.html"
    output:
        "{sample}/candidates/{sample}.all_candidate_variants.igv_report.html"
    conda: "envs/igvreport.yaml"
    threads: 1
    benchmark: "{sample}/benchmarks/combine_reports.tsv"
    shell: """ 
        python scripts/combine_igv_reports.py --sample {wildcards.sample} --snv {input.snv} --sv {input.sv} --cnv {input.cnv} --output {output}
    """
