# GAF annotation experiment for the HPRC 2023 Annual Meeting

This repo describes an experiment to use GAFs to represent annotations in a pangenome graph.

A genomic range can be represented as a path in the pangenome. 
The [Graph Alignment Format (GAF) text format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf), which was proposed to represent alignments, could be used to represent any type of annotation in a pangenome graph. 
To explore this approach within the vg toolkit, two subcommands were updated: 

- `vg gamsort` to sort and index bgzipped GAF files
- `vg annotate` to project annotation on the latest HPRC pangenomes.

Disclaimer: this is still experimental. The goal is to see what we can do already with minimal development. Multiple parts of the analysis below could be improved to make the whole process easier and more useful.

## Tools/dependencies

- [vg](https://github.com/vgteam/vg)
    - The *gafidx* branch contains code to sort/index GAF files and annotate BED/GFF3 on the HPRC pangenomes
    - A docker container with this version is available at: `quay.io/jmonlong/vg:gafidx`
- [sequenceTubeMap](https://github.com/vgteam/sequenceTubemap)
    - The *gafsupport* branch contains preliminary code to make sequenceTubeMap read GAF files.
    - A docker container with this version is available at: `quay.io/jmonlong/sequencetubemap:gaf`
- [Bandage](https://github.com/rrwick/Bandage)
    - We used v2022.09 of the forked version that can handle paths in GAF files: https://github.com/asl/BandageNG
- [Snakemake](snakemake.readthedocs.io/)
    - The pipeline was implemented in Snakemake. It uses the vg docker container mentioned above.
- Python3 to run the pre-processing script.

## Annotation data

For this experiment, we worked with annotations produced for each assembled haplotype of the HPRC freeze 1.
These are the haplotypes present in the pangenome (see below).
An index of the paths for each annotation file is available at [https://github.com/human-pangenomics/HPP_Year1_Assemblies/tree/main/annotation_index](https://github.com/human-pangenomics/HPP_Year1_Assemblies/tree/main/annotation_index).

We tested the pipeline on:

- BED files for the segmental duplication, tandem repeat, and RepeatMasker annotations
- GFF3 files for the gene annotation (from CAT)

```sh
wget https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/annotation_index/Year1_assemblies_v2_genbank_Seg_Dups.index
wget https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/annotation_index/Year1_assemblies_v2_genbank_TRF.index
wget https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/annotation_index/Year1_assemblies_v2_genbank_Repeat_Masker.index
wget https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/annotation_index/Year1_assemblies_v2_genbank_CAT_genes.index
```

## Pangenome

We used the CHM13-based freeze 1 Minigraph-Cactus pangenome.
It contains all the haplotypes as paths, has base-level resolution, and was known to be compatible with the vg toolkit.

The pangenome is available at [https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/freeze/freeze1/minigraph-cactus/](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes/freeze/freeze1/minigraph-cactus/).

The different indexes will be downloaded by the Snakemake pipeline.
To download them manually:

```sh
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-chm13.xg
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-chm13.gbwt
```

## Workflow

The workflow was implemented in Snakemake to annotate the different annotations on all the haplotypes. 
Under the hood, the workflow:

1. prepares the pangenomes
2. pre-process the annotation files
3. project the annotation onto the pangenome (`vg annotate`)
4. sort and index the pangenomic annotations (`vg gamsort`)

For example, assuming the full pangenome GBZ and raw annotations are ready, the commands might look like:

```sh 
## add the haplotype as a refernece path in the GBZ
vg gbwt -Z --set-tag "reference_samples=HG00438" --gbz-format -g gbz/HG00438.hprc-v1.0-mc-chm13.gbz hprc-v1.0-mc-chm13.gbz
## pre-process annotation
python3 pre-vg-annotate.py -i raw_annotation/HG00438.1.chm13.gff3.gz --add-prefix HG00438#1# --use-name-id | gzip > prep_annotation/gene_HG00438.1.gff3.gz
## project annotation onto pangenome
gunzip -c prep_annotation/gene_HG00438.1.gff3.gz | vg annotate -x gbz/HG00438.hprc-v1.0-mc-chm13.gbz -f - | vg convert -G - gbz/HG00438.hprc-v1.0-mc-chm13.gbz | gzip > unsorted_gaf/gene_HG00438.1.gaf.gz
## sort GAF
vg gamsort -t 1 -pG unsorted_gaf/gene_HG00438.1.gaf.gz | bgzip > gaf/gene_HG00438.1.gaf.gz
## index GAF
tabix -p gaf gaf/gene_HG00438.1.gaf.gz
```

To run the workflow:

```sh
snakemake -p --use-singularity
```

Of note, ressources to use SLURM are included in the workflow. 
Add `--slurm --latency-wait 30 --jobs 10`, for example, to the snakemake command above to run in SLURM systems.

## Pre-processing of the annotation files

The annotation files were modified using a python script ([`pre-vg-annotate.py`](pre-vg-annotate.py)) before running `vg annotate` to:

- Rename the elements so that it includes the haplotype of origin. This will later help differentiating the same genes/repeats/etc from different haplotypes.
- Create more descriptive names (for repeats).
- Add a prefix to the contig names to match the haplotype names in the pangenome (usually adding `{SAMPLE}#{HAP}#`).

### Pre-process gene annotation from CAT

```sh
python3 pre-vg-annotate.py -i {input} --add-prefix {pref} --use-name-id | gzip > {output}
```

where:

- `{pref}` is the haplotype name prefix (e.g. `HG00438#1#`) to add to the contig names.
- `--use-name-id` to replace the *Name* value with the value of *Name* and *ID*, seprated by a `:`. For example, the CAT annotations have `ID=HG00438.1_G0000001;Name=WASHC1` and the script will change it to `ID=HG00438.1_G0000001;Name=WASHC1:HG00438.1_G0000001`. 

### Pre-process tandem repeats from TRF

```sh
python3 pre-vg-annotate.py -i {input} --add-suffix "{suff}" --add-rep-n | gzip > {output}
```

where: 

-`{suff}` is the haplotype name suffix (e.g. `#HG00438#1`) to add to the repeat name.
- `--add-rep-n` to format the repeat names as `(<MOTIF>)<N>`

### Pre-process segmental duplications

```sh
python3 pre-vg-annotate.py -i {input} --add-suffix "{suff}" --add-len-fracm | gzip > {output}
```

where: 

-`{suff}` is the haplotype name suffix (e.g. `#HG00438#1`) to add to the SD name.
- `--add-len-fracm` to name the SD as `<LENGTH>bp_<FRACMATCH>`

### Pre-process RepeatMasker annotation

```sh
python3 pre-vg-annotate.py -i {input} --add-suffix "{suff}" --add-rm-class | gzip > {output}
```

where: 

-`{suff}` is the haplotype name suffix (e.g. `#HG00438#1`) to add to the SD name.
- `--add-rm-class` to prefix the repeat name with its class

## Visualization

### sequenceTubeMap

Using sequenceTubeMap, haplotypes, read alignments and paths can be visualized interactively. Hovering on a path displays its name, here the ID of a coding region of the BOLA2B gene. 

In the example below, the CHM13 haplotype is represented in purple, and the two haplotypes from HG00621 in greys.
The red and blue paths represents annotated coding sequence (CDS) for haplotype 1 and 2, respectively.

![](tubemap-cds.png)

### Bandage

Using BandageNG, a fork that can import paths in GAF files, paths were searched and colored to illustrate a mobile element insertion. 

In the example below, the CHM13 haplotype is colored in blue in purple, and an AluYa5 on haplotype 1 of HG00438 in yellow.
Thanks to the annotation we see that this non-reference sequence is likely a mobile element insertion.

![](bandage-mei.png)
