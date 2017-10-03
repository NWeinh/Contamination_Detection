# Detection of Cross-Sample Contamination in Next Generation Sequencing Data

### Analyze a single sample

For a single tumor/normal pair execute the contamination detection script:

```
Rscript cont.R -t <TUMOR BAM> -g <GERMLINE BAM> -n <SAMPLE ID> -p <PATH TO contPanel.csv> -o <OUTFILE>
```

### Analyze multiple samples (batch mode)

For multiple samples, modify file sampleSheet.csv and add sample IDs as well as location (full path) and name of tumor and control BAM files. One sample per row:

* SAMPLE_ID: 	patient identifier
* CONTROL: 	  BAM file name and location of the control (germline) sample
* TUMOR: 		  BAM file name and location of the tumor sample


```
Rscript cont.R -s <SAMPLESHEET> -p <contPanel.csv> -o <OUTFILE>
```

### List of available command line Options

```
        -t CHARACTER, --tumor=CHARACTER
                tumor bam file(s) [REQUIRED]

        -g CHARACTER, --germline=CHARACTER
                germline bam file(s) [REQUIRED]

        -p CHARACTER, --panel=CHARACTER
                panel of SNPs for tracking/contaminatin estimation [REQUIRED]

        -n CHARACTER, --sampleName=CHARACTER
                sample names(s) [REQUIRED]

        -s CHARACTER, --sampleSheet=CHARACTER
                sample sheet that lists sample name(s),
              tumor bam file location(s), germline bam file location(s)
              [REQUIRED]

        -e NUMERIC, --percentHom=NUMERIC
                percent homozygous [default= 10]

        -r NUMERIC, --minReads=NUMERIC
                minimum read depth at SNP position [default= 50]

        -c NUMERIC, --maxContLevelGerm=NUMERIC
                max contamination germline sample [default= 10]

        -q NUMERIC, --min_base_quality=NUMERIC
                minimum base quality [default= 20]

        -x NUMERIC, --contPerSNP=NUMERIC
                contamination per SNP [default= 0]

        -a NUMERIC, --aberrantSNP=NUMERIC
                number abberant SNPs [default= 5]

        -z NUMERIC, --aberrantSNPPercent=NUMERIC
                percent abberant SNPs [default= 10]

        -m CHARACTER, --mode=CHARACTER
                analysis mode [default= pair]

        -o CHARACTER, --out=CHARACTER
                output file name [default= conta.txt]

        -v LOGICAL, --verbose=LOGICAL
                get more detailed output [default= FALSE]

        -h, --help
                Show this help message and exit
```