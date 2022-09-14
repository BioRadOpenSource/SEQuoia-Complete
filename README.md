# SEQuoia-Complete

This Nextflow pipeline provides an easy to use analysis tool with
settings tailord for the Bio-Rad SEQuoia-Complete Kit.

### Preflight

You will need docker and nextflow installed:

```bash
# Install nextflow
wget -qO- https://get.nextflow.io | bash
```

The minimum recommended version is: 19.04.1

### Downloading Reference Genomes
It is suggested that you copy the tar of the reference that you want locally. These commands will take a while to run.
For full list of options see: [Sequoia Genomes](https://www.dropbox.com/sh/kqy6kt9qewqsmbl/AABSjlIs87-cWMLdLPd8eDOja?dl=0) 

When downloading from Dropbox, it will add `?dl=0` to the end of each link. This needs to be removed. This can be done manually, however using the output option in wget -O you can rename the file to match what is expected. In the use case of the hg38 genome, see the example below where the dropbox link with the `?dl=0` can be saved as the expected hg38.tar.gz with the -O option. 

** For the current version of complete only Human (hg38), Mouse (mm10), and Rat (rnor6) are supported **
 
```
mkdir ./ref_data/genome-annotations
cd ./ref_data/genome-annotations
wget -O hg38.tar.gz https://www.dropbox.com/s/hm6kyp70dtbqovr/hg38.tar.gz?dl=0
tar xvzf hg38.tar.gz
```

Then for your pipeline you will add the --genomes_base=./ref_data/genome-annotations

### Testing

Several tests come with this pipeline. To run time use the following
command. It requires thaty you have installed the reference genome as
above.

```bash
nextflow run . -w /mnt/nextflow/work/ -profile smoketest --genomes_base /mnt/ref_data/genome-annotations
```

It is recommended that you download the ref data as described above. If
not, data will be pulled from s3 as needed. 

### How to use this workflow

__Typical Usage:__

```bash
nextflow run . --reads './data/hg38/Sample11-273239772-sub/Sample11_Sub_S15_R1_001.fastq.gz'  --outDir /mnt/scratch/Sample11_Test --skipUmi -profile docker -w /mnt/nextflow/work/ --genomes_base /mnt/ref_data/genome-annotations --genome hg38 -resume --max_cpus 15 --max_memory 30
```

### Help

```bash
> nextflow run . --help
N E X T F L O W  ~  version 19.04.1
Launching `./main.nf` [cheeky_jennings] - revision: ae0e232204

Usage:

The typical command for running the pipeline is as follows:
nextflow run . --reads './tests/*_R{1,2}.fastq.gz' --genome hg38 --outDir /data/out --skipUmi --genomes_base /mnt/genome-annotations

Args:

REQUIRED:
    genome               (string )           Genome to align to and annotate against                                                                 [hg38, mm10, rnor6]                                      
    genomes_base         (string )           Bio-Rad formatted refence genomes and annotations                                                                                                                
    reads                (string )           Tese must be wrapped in single quotes. If R{1,2} is specified, UMI deduplication processes will be run.                     .*(R1_*\.fastq\.gz$|R2_*\.fastq\.gz$)

OPTIONAL:
    fivePrimeQualCutoff  (integer)           The read quality below which bases will be trimmed on the 5' end                                        [0, 42]                                                  
    max_cpus             (integer)           The max number of cpus the pipeline may use. Defaults provided by -profile.                                                                                      
    max_memory           (integer)           The max memory in GB that the pipeline may use. Defaults provided by -profile.                                                                                   
    minBp                (intger ) 15        Reads with fewer base pairs will be rejected                                                            [0, 500]                                                 
    minMapqToCount       (integer) 1         The minimum MapQ socre for an aligned read to count toward a feature count                              [0, 255]                                                 
    noTrim               (boolean)           Indicates whether or not trimming skipped on the reads                                                                                                           
    outDir               (string ) ./results Indicate the output directory to write to                                                                                                                        
    reverseStrand        (boolean)           Indicate if your library is reverse stranded                                                                                                                     
    skipUmi              (boolean)           Indicate that only R1 has been passed in and no UMI processing is required                                                                                       
    spikeType            (string ) NONE      The type of spike in used, if any                                                                       [NONE, ercc]                                             
    threePrimeQualCutoff (integer)           The read quality below which bases will be trimmed on the 3' end                                        [0, 42]                                                  
    validateInputs       (boolean) true      Ensure input meets standards and is below 500 million reads
```


## Contributing / Questions

Please feel free to make a PR or submit an issue. As indicated in the
licensing, all code is provided 'as is', but we will make an effort to
reply to issues and PR's when possible.
