<a href="https://rega.kuleuven.be/cev/viralmetagenomics"><img src="https://rega.kuleuven.be/cev/viralmetagenomics/pictures/lovm/image_preview" height="12.5%" width="12.5%" align="right"/></a>

# ViPER - Virome Paired-End Reads pipeline
[![Generic badge](https://img.shields.io/badge/GitHub-MatthijnssensLab-brightgreen?logo=github)](https://github.com/Matthijnssenslab)
[![Generic badge](https://img.shields.io/badge/NetoVIR-doi.org%2F10.1038%2Fsrep16532-blue)](https://doi.org/10.1038/srep16532)
[![Generic badge](https://img.shields.io/twitter/url?label=%40JMatthijnssens&style=social&url=https%3A%2F%2Ftwitter.com%2FJMatthijnssens)](https://twitter.com/JMatthijnssens)
[![Generic badge](https://img.shields.io/badge/Laboratory%20of%20Viral%20Metagenomics-1877F2?style=flat-square&logo=facebook&logoColor=white)](https://www.facebook.com/MatthijnssensLab)

The ViPER (Virome Paired-End Reads pipeline) script is used by the Laboratory of Viral Metagenomics to process raw paired-end Illumina reads. Reads are first trimmed by Trimmomatic, subsequently reads originating from the contaminome (sequenced and assembled negative controls) and the host genome can be removed by Bowtie2. Left-over reads are further assembled into contigs by metaSPAdes. To overcome the problem of viral genomes breaking into multiple pieces during assembly due to huge coverage, which makes the resulting De Bruijn graph too difficult to interpret by the assembler, a subset of 10 and 1% of the original reads may be applied (by `--triple-assembly`). These subsetted reads are also assembled and resulting contigs of all three assemblies (from original reads, 10% and 1% subset) are subsequently clustered together to remove redundancy in the contig set. This way, shorter contigs belonging to the same genome but from a different assembly will be removed and only the most complete contigs will be retained.

Final contigs can than be classified by DIAMOND and KronaTools with a lowest common ancestor approach. This is possible right after assembly with the `viper.sh` script, but is not necessary. Contigs can still be classified later on by the `viper-annotation.sh` script.

## Setup
```bash
git clone https://github.com/Matthijnssenslab/ViPER.git
cd ViPER
conda env create -f viper.yml
conda activate viper

#Install ViPER scripts in conda env bin
cd bin
ln -sfr * $CONDA_PREFIX/bin

#link Krona databases
rm -rf $CONDA_PREFIX/envs/viper/opt/krona/taxonomy
ln -sf /staging/leuven/stg_00029/DB/Krona/latest/ $CONDA_PREFIX/envs/viper/opt/krona/taxonomy
```

## Usage
```bash
usage: viper.sh [OPTIONS]

This script runs ViPER, the virome pipeline of the Laboratory of Viral Metagenomics (KU Leuven) for paired-end Illumina reads.

OPTIONS:
--------

REQUIRED:
   -1 | --read1			Path to the file with forward reads.
   -2 | --read2			Path to the file with reverse reads.
   
OPTIONAL:
 Trimming:
   -x | --crop			Crops reads with Trimmomatic CROP to this final length. First 19 bases of each read are removed by default with HEADCROP. (default:'')
   -p | --primer-file		Path to the primer file in fasta format with sequences that have to be trimmed by Trimmomatic, or a built-in option by Trimmomatic. 
   				(default: \$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa)

 Contamination removal:
   -c | --contaminome		Path to a bowtie2 indexed contaminome.
   -g | --host-genome		If specified, reads mapping to the given host genome will be removed. Requires the path to a bowtie2 indexed host genome.

 Assembly:
   -m | --min-length		The minimum length for final assembled scaffolds. (default: 500)
   -k | --spades-k-mer		List of k-mer sizes for SPAdes (must be odd and less than 128). (default: 21,33,55,77)
   --triple-assembly		Will perform three denovo assemblies with metaspades on the full reads, a 10% and 1% subset of the reads.
   				All assembled scaffolds will be concatenated and clustered together to remove redundancy (see also --cluster-cover/identity).
   --cluster-cover		% of the shortest sequence that should be covered during clustering. (default: 99)
   --cluster-identity		% of ANI for clustering scaffolds. (default: 99)
   --memory-limit		Memory to be reserved for SPAdes assembly in GB. (default: 250)

 Annotation:
   -d | --diamond-path		Path to diamond database. If not given, Diamond and KronaTools will be skipped.
   -s | --sensitivity		Can be 'default', 'fast', 'mid', 'more', 'very' and 'ultra' (default corresponds to --sensitive setting of Diamond).
   
GENERAL:
   -o | --outdir		Path where results will be stored and read files will be copied to (default: current directory). 
   -t | --threads		Number of threads to use. (default: 4)
   -h | --help    		Show this message and exit.
```


