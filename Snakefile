## read index files to match GENOME.[1/2] to each annotation file
def read_index(filen):
    inf = open(filen, 'rt')
    index_files = {}
    heads = next(inf).rstrip().split('\t')
    for line in inf:
        line = line.rstrip().split('\t')
        if 'reference' in heads and line[heads.index('reference')] != 'chm13':
            continue
        if line[heads.index('haplotype')] == '-':
            continue
        hap_name = '1' if line[heads.index('haplotype')] == 'paternal' else '2'
        hap_name = '{}.{}'.format(line[heads.index('sample')],
                                  hap_name)
        file_loc = line[heads.index('file_location')].replace('s3://human-pangenomics',
                                                             'https://s3-us-west-2.amazonaws.com/human-pangenomics')
        index_files[hap_name] = file_loc
    inf.close()
    return index_files

ann_files = {}
ann_files['gene'] = read_index('Year1_assemblies_v2_genbank_CAT_genes.index')
ann_files['trf'] = read_index('Year1_assemblies_v2_genbank_TRF.index')
ann_files['rm'] = read_index('Year1_assemblies_v2_genbank_Repeat_Masker.index')
ann_files['sd'] = read_index('Year1_assemblies_v2_genbank_Seg_Dups.index')

TYPES = ['gene', 'trf', 'sd', 'rm']
HAPS = list(ann_files['gene'].keys())

## run the workflow on just one haplotype
if 'one_hap' in config:
    HAPS = HAPS[0]

rule run:
    input:
        gaf=expand('gaf/{type}_{haplotype}.gaf.gz', type=TYPES, haplotype=HAPS),
        gaf_idx=expand('gaf/{type}_{haplotype}.gaf.gz.tbi', type=TYPES, haplotype=HAPS)

##
## TASKS
##

## sort and index the annotation GAFs
rule index_gaf:
    input: 'gaf/{type}_{haplotype}.gaf.gz'
    output: 'gaf/{type}_{haplotype}.gaf.gz.tbi'
    benchmark: 'benchmark/index_gaf.{type}_{haplotype}.benchmark.tsv'
    threads: 1
    resources:
        mem_mb="8GB",
        runtime="1h"
    singularity: 'docker://quay.io/jmonlong/vg:gafidx'
    shell: "tabix -p gaf {input}"

rule sort_gaf:
    input: 'unsorted_gaf/{type}_{haplotype}.gaf.gz'
    output: 'gaf/{type}_{haplotype}.gaf.gz'
    benchmark: 'benchmark/sort_gaf.{type}_{haplotype}.benchmark.tsv'
    threads: 1
    resources:
        mem_mb="8GB",
        runtime="3h"
    singularity: 'docker://quay.io/jmonlong/vg:gafidx'
    shell:
        """
        vg gamsort -t 1 -pG {input} | bgzip > {output}
        """

## make annotation GAFs
rule vg_annotate_bed:
    input:
        anno='prep_annotation/{anno_type}_{sample}.{haplotype}.bed.gz',
        gbz='gbz/{sample}.hprc-v1.0-mc-chm13.gbz'
    output: 'unsorted_gaf/{anno_type}_{sample}.{haplotype}.gaf.gz'
    log: 'log/vg_annotate.{anno_type}_{sample}.{haplotype}.log'
    benchmark: 'benchmark/vg_annotate.{anno_type}_{sample}.{haplotype}.benchmark.tsv'
    threads: 2
    resources:
        mem_mb="100GB",
        runtime="20h"
    singularity: 'docker://quay.io/jmonlong/vg:gafidx'
    shell: "gunzip -c {input.anno} | vg annotate -x {input.gbz} -b - 2> {log} | vg convert -G - {input.gbz} | gzip > {output}"
    
rule vg_annotate_gff3:
    input:
        anno='prep_annotation/{anno_type}_{sample}.{haplotype}.gff3.gz',
        gbz='gbz/{sample}.hprc-v1.0-mc-chm13.gbz'
    output: 'unsorted_gaf/{anno_type}_{sample}.{haplotype}.gaf.gz'
    log: 'log/vg_annotate.{anno_type}_{sample}.{haplotype}.log'
    benchmark: 'benchmark/vg_annotate.{anno_type}_{sample}.{haplotype}.benchmark.tsv'
    threads: 2
    resources:
        mem_mb="100GB",
        runtime="10h"
    singularity: 'docker://quay.io/jmonlong/vg:gafidx'
    shell: "gunzip -c {input.anno} | vg annotate -x {input.gbz} -f - 2> {log} | vg convert -G - {input.gbz} | gzip > {output}"    

## prepare the GBZ pangenome
rule prepare_haplotype_gbz:
    input: 'hprc{graph}.gbz'
    output: temp('gbz/{sample}.hprc{graph}.gbz')
    benchmark: 'benchmark/prep_haplotype_gbz.{sample}.{graph}.benchmark.tsv'
    threads: 1
    resources:
        mem_mb="25GB",
        runtime="1h"
    singularity: 'docker://quay.io/jmonlong/vg:gafidx'
    shell:
        """
        vg gbwt -Z --set-tag "reference_samples={wildcards.sample}" --gbz-format -g {output} {input}
        """

## prepare annotation files
rule prepare_cat_genes:
    input: 'raw_annotation/{haplotype}.chm13.gff3.gz'
    output: 'prep_annotation/gene_{haplotype}.gff3.gz'
    benchmark: 'benchmark/prep_cat_gene.{haplotype}.benchmark.tsv'
    threads: 1
    resources:
        mem_mb="8GB",
        runtime="3h"
    params:
        pref=lambda wildcards: wildcards.haplotype.replace('.', '#') + '#'
    shell: "python3 pre-vg-annotate.py -i {input} --add-prefix {params.pref} --use-name-id | gzip > {output}"

rule prepare_trf_repeats:
    input: 'raw_annotation/{haplotype}.trf.bed.gz'
    output: 'prep_annotation/trf_{haplotype}.bed.gz'
    benchmark: 'benchmark/prep_trf_repeats.{haplotype}.benchmark.tsv'
    threads: 1
    resources:
        mem_mb="8GB",
        runtime="3h"
    params:
        suff=lambda wildcards: "#" + wildcards.haplotype.replace('.', '#')
    shell:
        """
        python3 pre-vg-annotate.py -i {input} --add-suffix "{params.suff}" --add-rep-n | gzip > {output}
        """

rule prepare_segdups:
    input: 'raw_annotation/{haplotype}.sedef.bedpe'
    output: 'prep_annotation/sd_{haplotype}.bed.gz'
    benchmark: 'benchmark/prep_segdups.{haplotype}.benchmark.tsv'
    threads: 1
    resources:
        mem_mb="8GB",
        runtime="3h"
    params:
        suff=lambda wildcards: "#" + wildcards.haplotype.replace('.', '#')
    shell:
        """
        python3 pre-vg-annotate.py -i {input} --add-suffix "{params.suff}" --add-len-fracm | gzip > {output}
        """

rule prepare_repeat_masker:
    input: 'raw_annotation/{haplotype}.f1_assembly_v2_genbank_rm.bed'
    output: 'prep_annotation/rm_{haplotype}.bed.gz'
    benchmark: 'benchmark/prep_repeat_masker.{haplotype}.benchmark.tsv'
    threads: 1
    resources:
        mem_mb="8GB",
        runtime="3h"
    params:
        suff=lambda wildcards: "#" + wildcards.haplotype.replace('.', '#')
    shell:
        """
        python3 pre-vg-annotate.py -i {input} --add-suffix "{params.suff}" --add-rm-class | gzip > {output}
        """

## download annotation files
rule dwl_rm_annotation:
    output: 'raw_annotation/{haplotype}.f1_assembly_v2_genbank_rm.bed'
    threads: 1
    resources:
        mem_mb="4GB",
        runtime="1h"
    params:
        url=lambda wildcards: ann_files['rm'][wildcards.haplotype]
    shell: "wget -O {output} {params.url}"

rule dwl_cat_gene_annotation:
    output: 'raw_annotation/{haplotype}.chm13.gff3.gz'
    threads: 1
    resources:
        mem_mb="4GB",
        runtime="1h"
    params:
        url=lambda wildcards: ann_files['gene'][wildcards.haplotype]
    shell: "wget -O {output} {params.url}"

rule dwl_sd_annotation:
    output: 'raw_annotation/{haplotype}.sedef.bedpe'
    threads: 1
    resources:
        mem_mb="4GB",
        runtime="1h"
    params:
        url=lambda wildcards: ann_files['sd'][wildcards.haplotype]
    shell: "wget -O {output} {params.url}"

rule dwl_trf_annotation:
    output: 'raw_annotation/{haplotype}.trf.bed.gz'
    threads: 1
    resources:
        mem_mb="4GB",
        runtime="1h"
    params:
        url=lambda wildcards: ann_files['trf'][wildcards.haplotype]
    shell: "wget -O {output} {params.url}"

## download the pangenome
rule dwl_xg_gbwt:
    output:
        xg="hprc-v1.0-mc-chm13.xg",
        gbwt="hprc-v1.0-mc-chm13.gbwt"        
    threads: 1
    resources:
        mem_mb="4GB",
        runtime="1h"
    shell:
        """
        wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-chm13.xg
        wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-chm13.gbwt
        """

rule make_gbz:
    input:
        xg="hprc-v1.0-mc-chm13.xg",
        gbwt="hprc-v1.0-mc-chm13.gbwt"
    output: "hprc-v1.0-mc-chm13.gbz"
    params:
        paths_gbwt="temp_paths.gbwt",
        noref_gbwt="temp_noref.gbwt",
        comb_gbwt="temp_comb.gbwt"
    singularity: 'docker://quay.io/jmonlong/vg:gafidx'
    threads: 1
    resources:
        mem_mb="40GB",
        runtime="6h"
    shell:
        """
        vg gbwt -x {input.xg} -E -o {params.paths_gbwt}
        vg gbwt -R CHM13 -R _gbwt_ref -o {params.noref_gbwt} {input.gbwt}
        vg gbwt -m {params.noref_gbwt} {params.paths_gbwt} -o {params.comb_gbwt}
        vg gbwt -x {input.xg} {params.comb_gbwt} --gbz-format -g {output}
        rm -f {params.paths_gbwt} {params.noref_gbwt} {params.comb_gbwt}
        """
