#!/bin/bash
viper_command="$0 $@"
printf '%s\n' "$viper_command"

##### Setting default options #####
outdir="$PWD"
contaminome=''
spades_k_mer=21,33,55,77
diamond=0
diamond_path=''
threads=4
host_removal=0
host_genome=''
contaminome_removal=0
triple=0
trimmomatic_primer="$CONDA_PREFIX/share/trimmomatic/adapters/NexteraPE-PE.fa"
diamond_sensitivity="--sensitive"
minlength=500
cluster_cover=85
cluster_identity=95
crop=''
spades_memory=''
only_assembler=''
skip_trimming=0
move=1
read1_given=0
read2_given=0
unpaired=0
sample=''
warning=0

##### FUNCTIONS #####
#Help function
usage()
{
cat << EOF
usage: viper.sh [OPTIONS]

This script runs ViPER, the virome pipeline of the Laboratory of Viral Metagenomics (KU Leuven) for paired-end Illumina reads.

OPTIONS:
--------

REQUIRED:
   -1 | --read1			Path to the file with forward reads, may be gzipped.
   -2 | --read2			Path to the file with reverse reads, may be gzipped.
   
OPTIONAL:
 Reads:
   -u | --unpaired		File of unpaired reads if trimming was already performed beforehand (see --skip-trimming).
   
 Trimming:
   -x | --crop			Crops reads with Trimmomatic CROP to this final length. First 19 bases of each read are removed by default with HEADCROP. (default:'')
   -p | --primer-file		Path to the primer file in fasta format with sequences that have to be trimmed by Trimmomatic, or a built-in option by Trimmomatic. 
   				(default: \$CONDA_PREFIX/share/trimmomatic/adapters/NexteraPE-PE.fa)
   --skip-trimming		Continue with given reads and do not trim the reads for quality and adapters with Trimmomatic. 
   				Useful when you already have trimmed your reads beforehand with other software for example.

 Contamination removal:
   -c | --contaminome		Path to a bowtie2 indexed contaminome.
   -g | --host-genome		If specified, reads mapping to the given host genome will be removed. Requires the path to a bowtie2 indexed host genome.

 Assembly:
   -m | --min-length		The minimum length for final assembled contigs. (default: 500)
   -k | --spades-k-mer		List of k-mer sizes for SPAdes (must be odd and less than 128). (default: 21,33,55,77)
   --triple-assembly		Will perform three denovo assemblies with metaspades on the full reads, a 10% and 1% subset of the reads.
   				All assembled contigs will be concatenated and clustered together to remove redundancy (see also --cluster-cover/identity).
   --memory-limit		Memory (in GB) to be reserved for SPAdes assembly. (default: autodetected by SPAdes)
   --only-assembler		Runs only the assembler of metaspades, useful when the error correction gets stuck. (Caution: assembly without error correction does not yield the same results as normal assembly)

 Clustering (triple assembly enabled):
   --cluster-cover		% of the shortest sequence that should be covered during clustering. (default: 85)
   --cluster-identity		% of ANI for clustering contigs. (default: 95)
   --identify-proviruses	Automatically triggers --triple-assembly, identifies proviruses with genomad and performs contig integrity checks with CheckV. 
   				Proviruses are split into host and viral sequences, and contigs with problems identified by CheckV are separated before clustering.
				To reinclude these sequences you can run the separate script viper_cluster_study.py with the sequences from all your samples.
				See https://github.com/Matthijnssenslab/ViPER for more explanation.
   --checkv-db			Path to CheckV database. Required when --identify-proviruses is specified.
   --genomad-db			Path to genomad database. Required when --identify-proviruses is specified.

 Classification:
   -d | --diamond-db		Path to diamond database. If not given, Diamond and KronaTools will be skipped.
   -s | --sensitivity		Can be 'default', 'fast', 'mid', 'more', 'very' and 'ultra' (default corresponds to --sensitive setting of Diamond).
   
GENERAL:
   -h | --help    		Show this message.
   -o | --outdir		Path where results will be stored and read files will be copied to (default: current directory).
   -n | --name			Prefix to the output files, default is to use the common prefix of the read files or the date + timestamp, if no common prefix is found. 
   				Special characters are not allowed.
   -t | --threads		Number of threads to use. (default: 4)
   --bb-threads			Set the number of threads for BBMap tools (clumpify.sh, reformat.sh). (default: number of threads given to --threads)
				Can be set to one if viper.sh fails on these steps with 'Exception in thread' error.
   --keep-reads			Do not move the read files to the output directory, but keep them in place.

EOF
}

retfunc() {
return "$1"
}

#get file path
get_path() { 
echo "$(cd "$(dirname "$1")"; pwd -P)/$(basename "$1")"
}

#Get file names
get_name() {
file=$(basename -- "$1")
echo "$file"
}

#Find common prefix (and remove trailing dashes, dots, underscores and R)
common_prefix() {
printf '%s\n' "$1" "$2" | sed -e '1h;G;s,\(.*\).*\n\1.*,\1,;h;$!d' | sed -E -e 's/[-_\.]+R?[-_\.]*$//g'
}

check_error() {
	if [[ ! $? -eq 0 ]]; then
    	>&2 printf '\n%s\n' "[$(date "+%F %H:%M:%S")] ERROR: ViPER finished abnormally. Exit code: $1"
    	exit $1
	fi
}

#Check fasta file
#check_fasta() {
#perl -ne '
#    $id = />.+/;
#    die "Empty $.\n" if $id && $p || $id && eof;
#    $p = $id;
#    die "Invalid char $1 ($.)\n" if !$id && /([^A-Za-z\n])/
#    ' -- "$1"
#    }

check_fasta() {
seqkit stat "$1" | grep 'FASTA' &> /dev/null
}

##### CHECKS #####

#Check if all dependencies are installed in PATH
commands='seqkit samtools ktClassifyBLAST metaspades.py trimmomatic pigz bwa-mem2 diamond python bowtie2 reformat.sh fastqc perl clumpify.sh quast.py blastn makeblastdb anicalc.py aniclust.py'
for i in $commands; do
	command -v $i &> /dev/null
	if [[ ! $? -eq 0 ]]; then
    	printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: "$i" could not be found, please install "$i" in PATH or activate your conda environment."
    	exit 1
	fi
done

##### OPTIONS #####

if [[ $# -eq 0 ]]; then
	usage
	exit 1
fi

while [ ! $# -eq 0 ]; do
    case "$1" in
        -1 | --read1)
        	read1_given=1
        	if [[ -s "$2" ]]; then
        		read1_path=$(get_path "$2") 
        		shift
        	else
        		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The given forward read file does not exist or is empty."
        		exit 1
        	fi
        	;;
         -2 | --read2)
         	read2_given=1
        	if [[ -s "$2" ]]; then
        		read2_path=$(get_path "$2")
        		shift
        	else
        		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The given reverse read file does not exist or is empty."
        		exit 1
        	fi
        	;;
        -u | --unpaired)
        	if [[ -s "$2" ]]; then
        		unpaired_path=$(get_path "$2")
        		unpaired=1
        		shift
        	else
        		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The given unpaired read file does not exist or is empty."
        		exit 1
        	fi
        	;;
        -o | --outdir)
            if [[ "$2" != -* ]]; then
            	outdir=$(get_path "$2")
            	shift
            else
            	>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: You did not specify an output directory."
            	exit 1
            fi
            ;;
        -c | --contaminome)
            if [[ "$2" != -* ]]; then
            	contaminome_removal=1
                contaminome=$(get_path "$2")
                shift
            else
    			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: '-c | --contaminome' requires a bowtie2 indexed contaminome."
                exit 1
            fi
            ;;
        -m | --min-length)
        	if [[ "$2" =~ ^[0-9]+$ ]]; then
        		minlength=$2
        		shift
        	elif [[ "$2" == -* ]]; then
        		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Give a minimum length for assembled contigs."
        		exit 1
        	else
        		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The given minimum length is not an integer."
        		exit 1
        	fi
        	;;
        --memory-limit)
        	if [[ "$2" =~ ^[0-9]+$ ]]; then
        		spades_memory="-m $2"
        		shift
        	elif [[ "$2" == -* ]]; then
        		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Give a RAM limit for SPAdes assembly."
        		exit 1
        	else
        		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The given RAM limit is not an integer."
        		exit 1
        	fi
        	;;
		--only-assembler)
			only_assembler="--only-assembler"
        	;;
        -x | --crop)
        	if [[ "$2" =~ ^[0-9]+$ ]]; then
        		crop="CROP:$2"
        		shift
        	elif [[ "$2" == -* ]]; then
        		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Give a minimum length for assembled contigs."
        		exit 1
        	else
        		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The given minimum length is not an integer."
        		exit 1
        	fi
        	;;
        --skip-trimming)
        	skip_trimming=1
        	;;
        -g | --host-genome)
        	if [[ "$2" != -* ]]; then
        		host_removal=1
                host_genome=$(get_path "$2")
                shift
            else
                >&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: With '--host-genome' specified, viper.sh requires a bowtie2 indexed host genome."
                exit 1
            fi
            ;;
        -k | --spades-k-mer)
        	echo "$2" | \
				while read -d, i || [[ -n $i ]]; do 
					if [[ ! "$i" =~ ^[0-9]+$ ]]; then
						>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Given k-mers contain a non-integer."
						exit 1
					elif [ $((i%2)) -eq 0 ]; then 
						>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Even k-mer is not possible."
						exit 1
					elif [[ $i -gt 128 ]]; then
						>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: K-mer greater than 128."
						exit 1
					fi
				done
			if [[ $? -eq 0 ]]; then
        		spades_k_mer=$2
        		shift
        	else
        		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Something is wrong with k-mer list (see message above)."
        		exit 1
        	fi
        	;;
        -d | --diamond-db)
        	if [[ -s "$2" ]]; then
        		diamond=1
            	diamond_path=$(get_path "$2")
            	shift
            else
            	>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The provided diamond database does not exist."
            	exit 1
            fi
        	;;
        --triple-assembly)
        	triple=1
        	;;
        -p | --primer-file)
        	if [[ "$2" != -* ]]; then
        		check_fasta "$2"
        		if [[ $? -eq 0 ]]; then
                	trimmomatic_primer=$(get_path "$2")
                	shift
                else
                	>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The given primer file is not a valid fasta."
                	exit 1
                fi
            else
    			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: '--trimmomatic' requires a primer file in fasta format."
                exit 1
            fi
            ;;
        -s | --sensitivity)
        	if [[ "$2" != -* ]]; then
        		if [[ "$2" == "default" ]]; then
        			diamond_sensitivity="--sensitive"
        			shift
        		elif [[ "$2" == "mid" ]]; then
        			diamond_sensitivity="--mid-sensitive"
        			shift
        		elif [[ "$2" == "more" ]]; then
        			diamond_sensitivity="--more-sensitive"
        			shift
        		elif [[ "$2" == "very" ]]; then
        			diamond_sensitivity="--very-sensitive"
        			shift
        		elif [[ "$2" == "ultra" ]]; then
        			diamond_sensitivity="--ultra-sensitive"
        			shift
        		elif [[ "$2" == "fast" ]]; then
        			diamond_sensitivity=''
        			shift
        		else
        			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Unrecognized option for diamond sensitivity (options: default, fast, mid, more, very or ultra)."
        			exit 1
        		fi
        	else
        		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: --diamond-sensitivity requires an option (options: default, fast, mid, more, very or ultra)."
        		exit 1
            fi
            ;;
        --cluster-cover)
        	if [[ "$2" =~ ^[0-9]+$ ]]; then
        		cluster_cover=$2
        		shift
        	else
        		if [[ "$2" = -* ]]; then
        			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: No minimal % of covered sequence given for clustering."
        			exit 1
        		else
        			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Given value for --cluster-cover is not an integer."
        			exit 1
        		fi
        	fi
        	;;
        --cluster-identity)
        	if [[ "$2" =~ ^[0-9]+$ ]]; then
        		cluster_identity=$2
        		shift
        	else
        		if [[ "$2" = -* ]]; then
        			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: No % of sequence identity given for clustering."
        			exit 1
        		else
        			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Given value for --cluster-identity is not an integer."
        			exit 1
        		fi
        	fi
        	;;
		--identify-proviruses)
			identify_proviruses=1
			;;
		--checkv-db)
		    if [[ -s "$2" ]]; then
            	checkv_db=$(get_path "$2")
            	shift
            else
            	>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The provided CheckV database does not exist."
            	exit 1
            fi
			;;
		--genomad-db)
			if [[ -s "$2" ]]; then
            	genomad_db=$(get_path "$2")
            	shift
            else
            	>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The provided genomad database does not exist."
            	exit 1
            fi
        	;;
        -t | --threads)
        	if [[ "$2" =~ ^[0-9]+$ ]]; then
        		threads=$2
        		shift
        	else
        		if [[ "$2" = -* ]]; then
        			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] WARNING: No specified threads given, continuing with 4 threads."
					warning=$((warning+1))
        		else
        			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] WARNING: Given threads not an integer, continuing with 4 threads."
					warning=$((warning+1))
        			shift
        		fi
        	fi
        	;;
		--bb-threads)
        	if [[ "$2" =~ ^[0-9]+$ ]]; then
        		bbthreads=$2
        		shift
        	else
        		if [[ "$2" = -* ]]; then
					bbthreads=1
        			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] WARNING: No specified threads given, continuing with 1 thread for BBsuite tools."
					warning=$((warning+1))
        		else
					bbthreads=1
        			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] WARNING: Given threads not an integer, continuing with 1 thread for BBsuite tools."
					warning=$((warning+1))
        			shift
        		fi
        	fi
        	;;
        -n | --name)
        	if [[ "$2" == *['!'@#\$%^\&*()+]* ]]; then
        		:
        		shift
        	else
        		sample="$2"
        		shift
        	fi
        	;;
        --keep-reads)
        	move=0
        	;;
        -h | --help)
            usage
            exit
            ;;
        *)
            >&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: unrecognized option $1."
            usage
            exit 1
            ;;
    esac
    shift
done

#Set bbthreads
if [[ -z "$bbthreads" ]]; then
    bbthreads="$threads"
fi

#Check if required DBs are given when --identify-proviruses is set
if [[ $identify_proviruses -eq 1 ]]; then
	triple=1
	if [[ -d "$genomad_db" ]]; then
		if [[ ! -s  "$genomad_db/genomad_db" ]]; then
			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The database for genomad was not correctly provided."
			exit 1
		fi
	else
		if [[ ! -s "$genomad_db" ]]; then
			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The database for genomad was not correctly provided."
			exit 1
		fi
	fi

	if [[ ! -s "$checkv_db/hmm_db" ]]; then
		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The database for CheckV was not correctly provided."
		exit 1
	fi
fi


#Check if output directory already exists
if [[ -d "$outdir"/ASSEMBLY ]]; then
	>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The output directory already exists."
	exit 1
fi

### Check if all required options are given 
if [[ $read1_given -eq 1 ]]; then
	if [[ -s "$read1_path" ]]; then
		seqkit head "$read1_path" | seqkit stats | grep 'FASTQ' > /dev/null 2>&1
		if [[ ! $? -eq 0 ]]; then
			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The provided file "$read1_path" is not a FASTQ file."
			exit 1
		fi
	else
		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The provided path "$read1_path" does not lead to a file."
		exit 1
	fi
else
	>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: No forward reads given."
	exit 1
fi


if [[ $read2_given -eq 1 ]]; then
	if [[ -s "$read2_path" ]]; then
		seqkit head "$read2_path" | seqkit stats | grep 'FASTQ' > /dev/null 2>&1
		if [[ ! $? -eq 0 ]]; then
			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The provided file "$read2_path" is not a FASTQ file."
			exit 1
		fi
	else
		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The provided path "$read2_path" does not lead to a file."
		exit 1
	fi
else
	>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: No reverse reads given."
	exit 1
fi


if [[ $unpaired -eq 1 ]]; then
	if [[ $skip_trimming -eq 0 ]]; then
		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: You have to skip trimming if you want to provide an unpaired read file."
		exit 1
	else
		if [[ -s "$unpaired_path" ]]; then
			seqkit head "$unpaired_path" | seqkit stats | grep 'FASTQ' > /dev/null 2>&1
			if [[ ! $? -eq 0 ]]; then
				>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The provided file "$unpaired_path" is not a FASTQ file."
				exit 1
			fi
		else
			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The provided path "$unpaired_path" does not lead to a file."
			exit 1
		fi
	fi
fi

### Check if contaminome removal is specified and a valid indexed contaminome is given

if [[ $contaminome_removal -eq 1 ]]; then
	if [[ -n "$contaminome" ]]; then
		bowtie2-inspect -n "$contaminome" > /dev/null 2>&1
		if [[ ! $? -eq 0 ]]; then
			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The provided path by [-c | --contaminome] does not lead to a valid bowtie2 indexed contaminome."
			exit 1
		fi
	else
		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: No contaminome provided."
		exit 1
	fi
fi


### Check if given diamond database is valid 
if [[ $diamond -eq 1 ]]; then
	dbinfo=$(diamond dbinfo -p "$threads" --db "$diamond_path" --quiet | grep 'version' | grep -o -E [0-9]+) > /dev/null 2>&1
	if [[ ! $? -eq 0 ]]; then
		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The provided file is not a diamond database."
		exit 1
	elif [[ $dbinfo -le 1 ]]; then
		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: This database was made with an older version of diamond and is not compatible."
		>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Please remake your diamond database with a version of Diamond 2.0 or higher."
		exit 1
	fi
fi


### Check if host removal is specified and if there is a valid indexed genome provided

if [[ $host_removal -eq 1 ]]; then
	if [[ -z "$host_genome" ]]; then
		if [[ ! $? -eq 0 ]]; then
			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: No host genome specified."
			exit 1
		fi
	else
		bowtie2-inspect -n "$host_genome" > /dev/null 2>&1
		if [[ ! $? -eq 0 ]]; then
			>&2 printf '%s\n' "[$(date "+%F %H:%M:%S")] ERROR: The provided path does not lead to a valid bowtie2 indexed host genome."
			exit 1
		else 
			printf '%s\n' "[$(date "+%F %H:%M:%S")] INFO: Removing reads that map to $host_genome."
		fi
	fi
fi

# Extract file names
read1=$(get_name "$read1_path")
read2=$(get_name "$read2_path")
if [[ $unpaired -eq 1 ]]; then
	unpaired_name=$(get_name "$unpaired_path")
fi

if [[ -z "$sample" ]]; then
	sample=$(common_prefix "$read1" "$read2")
fi

#Test if there is a common prefix
if [[ -z "$sample" ]]; then
      >&2 printf '\n%s\n' "[$(date "+%F %H:%M:%S")] WARNING: No common prefix found between reads, continuing with date and timestamp as name. You might want to check if forward and reverse reads are from the same sample."
	  warning=$((warning+1))
      sample=$(date "+%Y%m%d-%H_%M_%S")
fi

##### START PIPELINE #####
printf '%s\n' "[$(date "+%F %H:%M:%S")] INFO: Starting ViPER!"

mkdir -p "$outdir"
cd "$outdir"

if [[ $move -eq 0 && $skip_trimming -eq 1 && $contaminome_removal -eq 0 && $host_removal -eq 0 && $triple -eq 0 ]]; then
	:
else
	mkdir -p READ
fi

if [[ $move -eq 1 ]]; then
	mv  "$read1_path" "$read2_path" READ/
	if [[ $unpaired -eq 1 ]]; then
		mv "$unpaired_path" READ/
		unpaired_path=$(get_path READ/"$unpaired_name")
	fi
	read1_path=$(get_path READ/"$read1")
	read2_path=$(get_path READ/"$read2")
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Moving reads to $(get_path READ)"
	cd READ
elif [[ $move -eq 0 && $skip_trimming -eq 1 && $contaminome_removal -eq 0 && $host_removal -eq 0 && $triple -eq 0 ]]; then
	:
else
	cd READ
fi

if [[ $skip_trimming -eq 0 ]]; then
	mkdir -p TRIMMED

### Trimming

	trimmomatic PE -threads "$threads" "$read1_path" "$read2_path" TRIMMED/"$sample".TRIM.R1.fastq.gz TRIMMED/"$sample".R1.unpaired.fastq.gz \
	TRIMMED/"$sample".TRIM.R2.fastq.gz TRIMMED/"$sample".R2.unpaired.fastq.gz \
	ILLUMINACLIP:"$trimmomatic_primer":2:30:7:1:true HEADCROP:19 LEADING:15 TRAILING:15 \
	SLIDINGWINDOW:4:20 MINLEN:50 $crop

	if [[ $? -eq 0 ]]; then
		printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Trimmomatic completed succesfully!"
	else
		>&2 printf '\n%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Something went wrong during trimming."
		exit 1
	fi

	cat TRIMMED/"$sample".R1.unpaired.fastq.gz TRIMMED/"$sample".R2.unpaired.fastq.gz > TRIMMED/"$sample".TRIM.unpaired.fastq.gz

	cd "$outdir"/READ/TRIMMED
	rm "$sample".R1.unpaired.fastq.gz
	rm "$sample".R2.unpaired.fastq.gz

	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Clumpifying trimmed reads for better compression."
	clumpify.sh reorder \
		in="$sample".TRIM.R1.fastq.gz \
		in2="$sample".TRIM.R2.fastq.gz \
		out="$sample".trimmed.R1.fastq.gz \
		out2="$sample".trimmed.R2.fastq.gz \
		ziplevel=9 \
		deleteinput=t \
		threads=$bbthreads
	clumpify.sh reorder \
		in="$sample".TRIM.unpaired.fastq.gz \
		out="$sample".trimmed.unpaired.fastq.gz \
		ziplevel=9 \
		deleteinput=t \
		threads=$bbthreads

	final_read1=$(get_path "$sample".trimmed.R1.fastq.gz)
	final_read2=$(get_path "$sample".trimmed.R2.fastq.gz)
	if [[ -s "$sample".trimmed.unpaired.fastq.gz ]]; then
		final_unpaired=$(get_path "$sample".trimmed.unpaired.fastq.gz)
		unpaired=1
	fi
else
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Skipped trimming as specified."
	final_read1="$read1_path"
	final_read2="$read2_path"
	if [[ $unpaired -eq 1 ]]; then
		final_unpaired="$unpaired_path"
	fi
fi
##############################################################################################################################################################

### Removing the contaminome 
#(you need to first make the de novo assembly of your negative controls and index it with bowtie2 to remove here)

if [[ $contaminome_removal -eq 1 ]]; then
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Removing contaminome."
	# Paired
	bowtie2 --very-sensitive -p "$threads" -x "$contaminome" \
	-1 "$final_read1" -2 "$final_read2" -S mapunmap_pair.sam
	if [[ ! $? -eq 0 ]]; then
		>&2 printf '\n%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Something went wrong during the contaminome removal of the paired reads."
		exit 1 
	fi
	samtools view -bS mapunmap_pair.sam -@ "$threads" | samtools view -@ "$threads" -b -f12 -F256 - | samtools sort -n - -o PEunmapped.sorted.bam -@ "$threads"
	samtools fastq PEunmapped.sorted.bam -1 NCout.R1.fastq -2 NCout.R2.fastq
	rm mapunmap_pair.sam
	rm PEunmapped.sorted.bam
	pigz -9 NCout.R1.fastq
	pigz -9 NCout.R2.fastq
	mv NCout.R1.fastq.gz "$sample".NCout.R1.fastq.gz
	mv NCout.R2.fastq.gz "$sample".NCout.R2.fastq.gz
	
	#Store names in variables
	final_read1=$(get_path "$sample".NCout.R1.fastq.gz)
	final_read2=$(get_path "$sample".NCout.R2.fastq.gz)

	if [[ $unpaired -eq 1 ]]; then
		# Unpaired
		bowtie2 --very-sensitive -p "$threads" -x "$contaminome" -U "$final_unpaired" -S mapunmap_unpair.sam
		if [[ ! $? -eq 0 ]]; then
			>&2 printf '\n%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Something went wrong during the contaminome removal of the unpaired reads."
			exit 1 
		fi
		samtools view -bS mapunmap_unpair.sam -@ "$threads" | samtools view -@ "$threads" -b -f4 -F256 - | samtools sort -n - -o UPunmapped.sorted.bam -@ "$threads"
		samtools fastq UPunmapped.sorted.bam > NCout.unpaired.fastq
		rm mapunmap_unpair.sam
		rm UPunmapped.sorted.bam 
		pigz -9 NCout.unpaired.fastq
		mv NCout.unpaired.fastq.gz "$sample".NCout.unpaired.fastq.gz
		
		#Store names in variables
		final_unpaired=$(get_path "$sample".NCout.unpaired.fastq.gz)
	fi
fi

##############################################################################################################################################################

### Removing host

if [[ $host_removal -eq 1 ]]; then
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Removing host genome."
	# Paired
	bowtie2 --very-sensitive -p "$threads" -x "$host_genome" -1 "$final_read1" -2 "$final_read2" -S mapunmap_pair.sam
	if [[ ! $? -eq 0 ]]; then
		>&2 printf '\n%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Something went wrong during the host genome removal of the paired reads."
		exit 1 
	fi
	samtools view -bS mapunmap_pair.sam -@ "$threads" | samtools view -@ "$threads" -b -f12 -F256 - | samtools sort -n - -o PEunmapped.sorted.bam -@ "$threads"
	samtools fastq PEunmapped.sorted.bam -1 Hostout.R1.fastq -2 Hostout.R2.fastq -@ "$threads"
	rm mapunmap_pair.sam
	rm PEunmapped.sorted.bam
	pigz -9 Hostout.R1.fastq
	pigz -9 Hostout.R2.fastq
	mv Hostout.R1.fastq.gz "$sample".Hostout.R1.fastq.gz
	mv Hostout.R2.fastq.gz "$sample".Hostout.R2.fastq.gz
	
	#Store names in variables
	final_read1=$(get_path "$sample".Hostout.R1.fastq.gz)
	final_read2=$(get_path "$sample".Hostout.R2.fastq.gz)

	if [[ $unpaired -eq 1 ]]; then
		# Unpaired
		bowtie2 --very-sensitive -p "$threads" -x "$host_genome" -U "$final_unpaired" -S mapunmap_unpair.sam
		if [[ ! $? -eq 0 ]]; then
			>&2 printf '\n%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Something went wrong during the host genome removal of the unpaired reads."
			exit 1 
		fi
		samtools view -bS mapunmap_unpair.sam -@ "$threads" | samtools view -@ "$threads" -b -f4 -F256 - | samtools sort -n - -o UPunmapped.sorted.bam -@ "$threads"
		samtools fastq UPunmapped.sorted.bam > Hostout.unpaired.fastq -@ "$threads"
		rm mapunmap_unpair.sam
		rm UPunmapped.sorted.bam
		pigz -9 Hostout.unpaired.fastq
		mv Hostout.unpaired.fastq.gz "$sample".Hostout.unpaired.fastq.gz
		
		#Store names in variables
		final_unpaired=$(get_path "$sample".Hostout.unpaired.fastq.gz)
	fi
fi

##############################################################################################################################################################

### QC of reads
printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Checking read quality with fastqc."
mkdir -p "$outdir"/QC/FASTQC
if [[ $unpaired -eq 1 ]]; then
	fastqc -o "$outdir"/QC/FASTQC -t "$threads" -q "$final_read1" "$final_read2" "$final_unpaired"
else
	fastqc -o "$outdir"/QC/FASTQC -t "$threads" -q "$final_read1" "$final_read2"
fi

# After all samples are done you can run multiqc to output the QC of all samples in 1 file
#multiqc -o QC .

##############################################################################################################################################################

### Subsetting the clean reads with reformat.sh from the BBmap suite

if [[ $triple -eq 1 ]]; then
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Subsetting reads for triple assembly."
	# 10% of reads
	reformat.sh in="$final_read1" in2="$final_read2" \
		out1="$sample".subset_10.R1.fastq.gz out2="$sample".subset_10.R2.fastq.gz \
		samplerate=0.1 sampleseed=1234 threads=$bbthreads

	# 1% of reads
	reformat.sh in="$final_read1" in2="$final_read2" \
		out1="$sample".subset_1.R1.fastq.gz out2="$sample".subset_1.R2.fastq.gz \
		samplerate=0.01 sampleseed=1234 threads=$bbthreads

	subset10_R1=$(get_path "$sample".subset_10.R1.fastq.gz)
	subset10_R2=$(get_path "$sample".subset_10.R2.fastq.gz)
	subset1_R1=$(get_path "$sample".subset_1.R1.fastq.gz)
	subset1_R2=$(get_path "$sample".subset_1.R2.fastq.gz)
fi
	
##############################################################################################################################################################

### Assemblies
if [[ $triple -eq 1 ]]; then
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Starting triple assembly with metaSPAdes."
	# Full assembly
	mkdir -p "$outdir"/ASSEMBLY
	cd "$outdir"/ASSEMBLY
	
	if [[ $unpaired -eq 1 ]]; then
		metaspades.py -1 "$final_read1" -2 "$final_read2" \
		-s "$final_unpaired" -t "$threads" -k "$spades_k_mer" -o ASSEMBLY1 $spades_memory $only_assembler
	else
		metaspades.py -1 "$final_read1" -2 "$final_read2" \
		-t "$threads" -k "$spades_k_mer" -o ASSEMBLY1 $spades_memory $only_assembler
	fi

	check_error 2
	
	cd ASSEMBLY1
	mv contigs.fasta "$sample".full.contigs.fasta
	#to add the sample names to your assemblies
	sed -i "s/NODE_/NODE_A/g" "$sample".full.contigs.fasta
	sed -i "s/>.*/&_${sample}/" "$sample".full.contigs.fasta

	# Rename scaffolds too
	mv scaffolds.fasta "$sample".full.scaffolds.fasta
	#to add the sample names to your assemblies
	sed -i "s/NODE_/NODE_A/g" "$sample".full.scaffolds.fasta
	sed -i "s/>.*/&_${sample}/" "$sample".full.scaffolds.fasta

	# 10% assembly
	printf '\n\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Starting second assembly with 10% of the reads."
	cd "$outdir"/ASSEMBLY
	metaspades.py -1 "$subset10_R1" -2 "$subset10_R2" \
	-t "$threads" -k "$spades_k_mer" -o ASSEMBLY2 $spades_memory
	
	check_error 3

	cd ASSEMBLY2
	mv contigs.fasta $"$sample".10-percent.contigs.fasta
	#to add the sample names to your assemblies
	sed -i "s/NODE_/NODE_B/g" "$sample".10-percent.contigs.fasta
	sed -i "s/>.*/&_${sample}/" "$sample".10-percent.contigs.fasta

	# Rename scaffolds too
	mv scaffolds.fasta $"$sample".10-percent.scaffolds.fasta
	#to add the sample names to your assemblies
	sed -i "s/NODE_/NODE_B/g" "$sample".10-percent.scaffolds.fasta
	sed -i "s/>.*/&_${sample}/" "$sample".10-percent.scaffolds.fasta

	# 1% assembly
	printf '\n\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Starting third assembly with 1% of the reads."
	cd "$outdir"/ASSEMBLY
	metaspades.py -1 "$subset1_R1" -2 "$subset1_R2" \
	-t "$threads" -k "$spades_k_mer" -o ASSEMBLY3 $spades_memory
	
	check_error 4
	
	cd ASSEMBLY3
	mv contigs.fasta "$sample".1-percent.contigs.fasta
	#to add the sample names to your assemblies
	sed -i "s/NODE_/NODE_C/g" "$sample".1-percent.contigs.fasta
	sed -i "s/>.*/&_${sample}/" "$sample".1-percent.contigs.fasta

	# Rename scaffolds too
	mv scaffolds.fasta "$sample".1-percent.scaffolds.fasta
	#to add the sample names to your assemblies
	sed -i "s/NODE_/NODE_C/g" "$sample".1-percent.scaffolds.fasta
	sed -i "s/>.*/&_${sample}/" "$sample".1-percent.scaffolds.fasta

else
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Starting assembly with metaSPAdes."
	cd "$outdir"
	if [[ $unpaired -eq 1 ]]; then
		metaspades.py -1 "$final_read1" -2 "$final_read2" \
		-s "$final_unpaired" -t "$threads" -k "$spades_k_mer" -o ASSEMBLY $spades_memory
	else
		metaspades.py -1 "$final_read1" -2 "$final_read2" \
		-t "$threads" -k "$spades_k_mer" -o ASSEMBLY $spades_memory
	fi
	cd ASSEMBLY
	mv contigs.fasta "$sample".contigs.fasta
	#to add the sample names to your assemblies
	sed -i "s/>.*/&_${sample}/" "$sample".contigs.fasta

	# Rename scaffolds too
	mv scaffolds.fasta "$sample".scaffolds.fasta
	#to add the sample names to your assemblies
	sed -i "s/>.*/&_${sample}/" "$sample".scaffolds.fasta
fi
printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Assembly finished!"
##############################################################################################################################################################

### Clustering assemblies
cd "$outdir"
mkdir -p CONTIGS
cd CONTIGS

if [[ $triple -eq 1 ]]; then
	mkdir -p triple-assembly
	cp "$outdir"/ASSEMBLY/ASSEMBLY*/*.contigs.fasta triple-assembly/
	cat triple-assembly/*.contigs.fasta > "$sample"_all.contigs.fasta
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Filtering contigs larger than "$minlength"bp."
	seqkit seq -m "$minlength" -j "$threads" "$sample"_all.contigs.fasta > "$sample"_"$minlength"-unclustered.contigs.fasta
	seqkit sort --by-length --reverse -o "$sample"_"$minlength"-unclustered.contigs.fasta "$sample"_"$minlength"-unclustered.contigs.fasta
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Clustering contigs on "$cluster_identity"% identity over "$cluster_cover"% of the length."
	#Cluster_genomes.pl -f "$sample"_"$minlength"-unclustered.contigs.fasta -c "$cluster_cover" -i "$cluster_identity" -t "$threads"

	if [[ $identify_proviruses -eq 1 ]]; then
		viper_cluster.py -i "$sample"_"$minlength"-unclustered.contigs.fasta -o "$sample"_"$minlength" -t $threads --checkv-db "$checkv_db" \
									--genomad-db "$genomad_db" --min-identity $cluster_identity --min-coverage $cluster_cover
	else
		viper_cluster.py -i "$sample"_"$minlength"-unclustered.contigs.fasta -o "$sample"_"$minlength" -t $threads \
									--min-identity $cluster_identity --min-coverage $cluster_cover
	fi

	mv "$sample"_"$minlength".fasta "$sample"_"$minlength".contigs.fasta

else
	cp "$outdir"/ASSEMBLY/"$sample".contigs.fasta .
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Filtering contigs larger than "$minlength"bp."
	seqkit seq -m "$minlength" -j "$threads" "$sample".contigs.fasta > "$sample"_"$minlength".contigs.fasta
fi

contigs="$sample"_"$minlength".contigs.fasta

printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Running Quast for assembly QC."
quast.py -t "$threads" --fast -m "$minlength" -o "$outdir"/QC/QUAST "$outdir"/CONTIGS/"$contigs"

if [[ ! $? -eq 0 ]]; then
	quast=1
else
	quast=0
fi
##############################################################################################################################################################

### Diamond taxonomical annotation
if [[ $diamond -eq 1 ]]; then
	cd "$outdir"
	mkdir -p DIAMOND
	cd DIAMOND
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Running Diamond!"
	diamond blastx --db "$diamond_path" -q "$outdir"/CONTIGS/"$contigs" -a "$sample" -p "$threads" $diamond_sensitivity -c 1 -b 5 --tmpdir /dev/shm
	diamond view -a "$sample" -o "$sample".m8 -p "$threads"
	
	if [[ ! $? -eq 0 ]]; then
		>&2 printf '\n%s\n' "[$(date "+%F %H:%M:%S")] ERROR: Something went wrong with Diamond."
		exit 1 
	fi
fi
############################################################################################################################################################

if [[ $diamond -eq 1 ]]; then
### Relative abundances by mapping 
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Counting abundances for Krona."
	cd "$outdir"/CONTIGS

	bwa-mem2 index "$contigs"
	bwa-mem2 mem "$contigs" "$final_read1" "$final_read2" -t "$threads" | samtools view -Su - | samtools sort - -o "$sample".R.sort.bam

	if [[ $unpaired -eq 1 ]]; then
		bwa-mem2 mem "$contigs" "$final_unpaired" -t "$threads" | samtools view -Su - | samtools sort - -o "$sample".un.sort.bam
		samtools merge -f "$sample".bam "$sample".R.sort.bam "$sample".un.sort.bam
		rm "$sample".R.sort.bam
		rm "$sample".un.sort.bam
	else
		mv "$sample".R.sort.bam "$sample".bam
	fi
	
	samtools index "$sample".bam
	samtools idxstats "$sample".bam | cut -f1,3 > "$sample".magnitudes

### Krona visualization
	cd "$outdir"
	mkdir -p KRONA
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: Making Krona chart."
	ktImportBLAST -o KRONA/"$sample".html "$outdir"/DIAMOND/"$sample".m8,"$sample" "$outdir"/DIAMOND/"$sample".m8:"$outdir"/CONTIGS/"$sample".magnitudes,"$sample".magn

	#ktClassifyBLAST "$outdir"/DIAMOND/"$sample".m8 -o KRONA/"$sample".tab
	#awk 'NR==FNR { a[$1]=$2; next} $1 in a {print $0,"\t"a[$1]}' "$outdir"/CONTIGS/"$sample".magnitudes "$outdir"/KRONA/"$sample".tab > "$outdir"/KRONA/"$sample".magnitudes.tab
elif [[ ! $quast -eq 0 ]]; then
	retfunc 1
else
	retfunc 0
fi

if [[ ! $? -eq 0 ]]; then
	check_error 5
elif [[ $warning -gt 0 ]]; then
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: ViPER finished with $warning warning(s)."
else
	printf '\n%s\n' "[$(date "+%F %H:%M:%S")] INFO: ViPER finished successfully!"
fi
