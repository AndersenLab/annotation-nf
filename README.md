# annotation-nf

Annotate VCF with snpeff and bcsq


# Pipeline overview

```
-------------    
ANNOTATION-NF
-------------

nextflow main.nf --debug

nextflow main.nf --vcf=hard-filtered.vcf --species=c_elegans --divergent_regions=divergent_regions_strain.bed

    parameters           description                                              Set/Default
    ==========           ===========                                              ========================
    --debug              Set to 'true' to test                                    ${params.debug}
    --species            Species: 'c_elegans', 'c_tropicalis' or 'c_briggsae'     ${params.species}
    --vcf                hard filtered vcf to calculate variant density           ${params.vcf}
    --divergent_regions  (Optional) Divergent region bed file                     ${params.divergent_regions}
    --reference          Reference used based on species and project              ${params.reference}
    --output             (Optional) output folder name                            ${params.output}
 
    username                                                                      ${"whoami".execute().in.text}

    HELP: http://andersenlab.org/dry-guide/pipeline-annotation-nf
```

![Pipeline-overview](bin/annotation-nf.drawio.svg)

## Software Requirements

* The latest update requires Nextflow version 20.0+. On QUEST, you can access this version by loading the `nf20_env` conda environment prior to running the pipeline command:

```
module load python/anaconda3.6
source activate /projects/b1059/software/conda_envs/nf20_env
```

Alternatively you can update Nextflow by running:

```
nextflow self-update
```


# Usage

## Testing on Quest

*This command uses a test dataset*

```
nextflow run andersenlab/annotation-nf --debug
```

## Running on Quest

You should run this in a screen session.

```
nextflow run andersenlab/annotation-nf --vcf <path_to_vcf> --species <species> --divergent_regions <path_to_file>
```

# Parameters

## --debug

You should use `--debug true` for testing/debugging purposes. This will run the debug test set (located in the `test_data` folder).

For example:

```
nextflow run andersenlab/annotation-nf --debug
```

## --vcf

Path to the hard-filter, isotype VCF (output from `post-gatk-nf`)

## --species

Choose from c_elegans, c_briggsae, or c_tropicalis. Species will specifiy a default reference genome. You can select a different one if you prefer (see below)

## --divergent_regions

This is the `divergent_regions_strain.bed` file output from the `post-gatk-nf` pipeline. This file is used to add a column to the flat file if the variant is within a divergent region. Currently, *C. elegans* is the only species with divergent regions, if running for another species, do not provide a divergent_regions file and the pipeline will ignore it.

### --reference, --project, --ws_build (optional)

By default, the reference genome is set by the species parameter. If you don't want to use the default, you could change the project and/or ws_build. As long as the genome is in the proper location on quest (for more, see the [genomes-nf](pipeline-genomes-nf.md) pipeline), this will work. Alternatively, you could provide the path to a reference of your choice.

**Defaults:**
- *C. elegans* - `/projects/b1059/data/c_elegans/genomes/PRJNA13758/WS276/c_elegans.PRJNA13758.WS276.genome.fa.gz`
- *C. briggsae* - `/projects/b1059/data/c_briggsae/genomes/QX1410_nanopore/Feb2020/c_briggsae.QX1410_nanopore.Feb2020.genome.fa.gz`
- *C. tropicalis* - `/projects/b1059/data/c_tropicalis/genomes/NIC58_nanopore/June2021/c_tropicalis.NIC58_nanopore.June2021.genome.fa.gz`

### --ncsq_param (optional)

This parameter is necessary for correct annotation using BCSQ for variants with many different annotations (like found in divergent regions). In 20210121 we found that the default value of `224` was sufficient, but as more strains are added this number might need to increase. If there is an issue, you should see a warning error from BCFtools and they should suggest what to change this parameter to.

# Output

```
├── strain_vcf
│   ├── {strain}.{date}.vcf.gz
│   └── {strain}.{date}.vcf.gz.tbi
└── variation
    ├── WI.{date}.hard-filter.isotype.snpeff.vcf.gz
    ├── WI.{date}.hard-filter.isotype.snpeff.vcf.gz.tbi
    ├── snpeff.stats.csv
    ├── WI.{date}.hard-filter.isotype.bcsq.vcf.gz
    ├── WI.{date}.hard-filter.isotype.bcsq.vcf.gz.tbi
    ├── WI.{date}.hard-filter.isotype.bcsq.vcf.gz.stats.txt 
    └── WI.{date}.strain-annotation.bcsq.tsv
 
```

