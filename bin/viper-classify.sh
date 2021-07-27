#!/bin/bash

##### Setting default options #####
outdir="$PWD"
diamond=0
diamond_path=''
threads=4
diamond_sensitivity="--sensitive"
minlength=500
filterlength=0

##### FUNCTIONS #####
#Help function
usage()
{
cat << EOF
usage: viper-classify.sh [OPTIONS]

This script takes a fasta file as input and uses Diamond and Krona to classify the input sequences and make a Krona pie chart.

OPTIONS:
--------

REQUIRED:
   -f | --fasta			Input fasta file with contigs to annotate.
   -1 | --read1			Path to the fastq file with forward reads.
   -2 | --read2			Path to the fastq file with reverse reads.
   -d | --diamond-path		Path to diamond database.
   
OPTIONAL:
   -u | --unpaired		Path to the fastq file with unpaired reads.
   -m | --min-length		The minimum length for final assembled scaffolds. (default: 500)
   -s | --sensitivity		Can be 'default', 'fast', 'mid', 'more', 'very' and 'ultra' (default corresponds to --sensitive setting of Diamond).
   
GENERAL:
   -o | --outdir		Output directory (default: current directory). 
   -t | --threads		Number of threads to use. (default: 4)
   -h | --help    		Show this message and exit.

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
printf '%s\n' "$1" "$2" | sed -e '1h;G;s,\(.*\).*\n\1.*,\1,;h;$!d' | sed -E -e 's/[-_\.]+R?$//g'
}

#Check fasta file
check_fasta() {
perl -ne '
    $id = />.+/;
    die "Empty $.\n" if $id && $p || $id && eof;
    $p = $id;
    die "Invalid char $1 ($.)\n" if !$id && /([^A-Za-z\n])/
    ' -- "$1"
    }

##### OPTIONS #####

if [[ $# -eq 0 ]]; then
	usage
	exit 1
fi

while [ ! $# -eq 0 ]; do
    case "$1" in
    	-f | --fasta)
    		if [[ -s "$2" ]]; then
    			check_fasta "$2"
    			if [[ $? -eq 0 ]]; then
        			fasta=$(get_path "$2")
        			shift
        		else
        			>&2 printf '\n%s\n\n' "[ERROR]: The given input file is not a fasta file."
        			exit 1
        		fi
        	else
        		>&2 printf '\n%s\n\n' "[ERROR]: The given input file does not exist or is empty."
        		exit 1
        	fi
        	;;
        -1 | --read1)
        	if [[ -s "$2" ]]; then
        		read1_path=$(get_path "$2") 
        		shift
        	else
        		>&2 printf '\n%s\n\n' "[ERROR]: The given forward read file does not exist or is empty."
        		exit 1
        	fi
        	;;
         -2 | --read2)
        	if [[ -s "$2" ]]; then
        		read2_path=$(get_path "$2")
        		shift
        	else
        		>&2 printf '\n%s\n\n' "[ERROR]: The given reverse read file does not exist or is empty."
        		exit 1
        	fi
        	;;
        -u | --unpaired)
        	if [[ -s "$2" ]]; then
        		unpaired_path=$(get_path "$2") 
        		shift
        	else
        		>&2 printf '\n%s\n\n' "[ERROR]: The given unpaired read file does not exist or is empty."
        		exit 1
        	fi
        	;;
        -o | --outdir)
            if [[ "$2" != -* ]]; then
            	outdir=$(get_path "$2")
            	shift
            else
            	>&2 printf '\n%s\n\n' "[ERROR]: You did not specify an output directory."
            	exit 1
            fi
            ;;
        -m | --min-length)
        	if [[ "$2" =~ ^[0-9]+$ ]]; then
        		minlength=$2
        		filterlength=1
        		shift
        	elif [[ "$2" == -* ]]; then
        		>&2 printf '\n%s\n\n' "[ERROR]: Give a minimum length for assembled scaffolds."
        		exit 1
        	else
        		>&2 printf '\n%s\n\n' "[ERROR]: The given minimum length is not an integer."
        		exit 1
        	fi
        	;;
        -d | --diamond-path)
        	if [[ -s "$2" ]]; then
        		diamond=1
            	diamond_path=$(get_path "$2")
            	shift
            else
            	>&2 printf '\n%s\n\n' "[ERROR]: The provided diamond database does not exist."
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
        			>&2 printf '\n%s\n\n' "[ERROR]: Unrecognized option for diamond sensitivity (options: default, fast, mid, more, very or ultra)."
        			exit 1
        		fi
        	else
        		>&2 printf '\n%s\n\n' "[ERROR]: --diamond-sensitivity requires an option (options: default, fast, mid, more, very or ultra)."
        		exit 1
            fi
            ;;
        -t | --threads)
        	if [[ "$2" =~ ^[0-9]+$ ]]; then
        		threads=$2
        		shift
        	else
        		if [[ "$2" = -* ]]; then
        			>&2 printf '\n%s\n\n' "[WARNING]: No specified threads given, continuing with 4 threads."
        		else
        			>&2 printf '\n%s\n\n' "[WARNING]: Given threads not an integer, continuing with 4 threads."
        			shift
        		fi
        	fi
        	;;
        -h | --help)
            usage
            exit
            ;;
        *)
            >&2 printf '\n%s\n\n' "[ERROR]: unrecognized option $1."
            usage
            exit 1
            ;;
    esac
    shift
done

#Check if all dependencies are installed in PATH
commands='seqkit samtools ktClassifyBLAST bwa-mem2 diamond'
for i in $commands; do
	command -v $i &> /dev/null
	if [[ ! $? -eq 0 ]]; then
    printf '\n%s\n' "[ERROR]: "$i" could not be found, please install "$i" in PATH or activate your conda environment."
    exit 1
	fi
done

### Check if all required options are given 

if [[ -s "$read1_path" ]]; then
	seqkit head "$read1_path" | seqkit stats | grep 'FASTQ' > /dev/null 2>&1
	if [[ ! $? -eq 0 ]]; then
		>&2 printf '\n%s\n\n' "[ERROR]: The provided file "$read1_path" is not a FASTQ file."
		exit 1
	fi
else
	>&2 printf '\n%s\n\n' "[ERROR]: The provided path "$read1_path" does not lead to a file."
fi

if [[ -s "$read2_path" ]]; then
	seqkit head "$read2_path" | seqkit stats | grep 'FASTQ' > /dev/null 2>&1
	if [[ ! $? -eq 0 ]]; then
		>&2 printf '\n%s\n\n' "[ERROR]: The provided file "$read2_path" is not a FASTQ file."
		exit 1
	fi
else
	>&2 printf '\n%s\n\n' "[ERROR]: The provided path "$read2_path" does not lead to a file."
fi

### Check if given diamond database is valid 
if [[ $diamond -eq 1 ]]; then
	dbinfo=$(diamond dbinfo -p "$threads" -d "$diamond_path" --quiet | grep 'version' | grep -o -E [0-9]+) > /dev/null 2>&1
	if [[ ! $? -eq 0 ]]; then
		>&2 printf '\n%s\n\n' "[ERROR]: The provided file is not a diamond database."
		exit 1
	elif [[ $dbinfo -le 1 ]]; then
		>&2 printf '\n%s' "[ERROR]: This database was made with an older version of diamond and is not compatible."
		>&2 printf '\n%s\n\n' "[ERROR]: Please remake your diamond database with a version of Diamond 2.0 or higher."
		exit 1
	fi
else
	>&2 printf '\n%s\n\n' "[ERROR]: No diamond database given."
	exit 1
fi

# Extract file names
read1=$(get_name "$read1_path")
read2=$(get_name "$read2_path")
#sample=$(common_prefix "$read1" "$read2")

fasta_name=$(basename -- "$fasta")
sample="${fasta_name%%.*}"

#Test if there is a common prefix
#if [[ -z "$sample" ]]; then
#      >&2 printf '\n%s\n\n' "[WARNING]: No common prefix found between reads, continuing with generic sample name but you might want to check if forward and reverse reads are from the same sample."
#      sample="sample"
#fi

##### START ANNOTATION PIPELINE #####
printf '\n%s\n\n' "[INFO]: Starting ViPER annotation script!"
mkdir -p "$outdir"
cd "$outdir"


if [[ $filterlength -eq 1 ]]; then
	printf '\n%s\n\n' "[INFO]: Filtering sequences larger than "$minlength"bp."
	seqkit seq -m "$minlength" -j "$threads" "$fasta" > "$sample"_"$minlength".fa
	seqkit sort --by-length --reverse -o "$sample"_"$minlength".fa "$sample"_"$minlength".fa
	fasta="$sample"_"$minlength".fa
fi
	
	
### Diamond taxonomical annotation
if [[ $diamond -eq 1 ]]; then
	printf '\n%s\n\n' "[INFO]: Running Diamond!"
	diamond blastx -d "$diamond_path" -q "$fasta" -a "$sample" -p "$threads" $diamond_sensitivity -c 1 -b 5 --tmpdir /dev/shm
	diamond view -d "$diamond_path" -a "$sample" -o "$sample".m8 -p "$threads"
	
	if [[ ! $? -eq 0 ]]; then
		>&2 printf '\n%s\n\n' "[ERROR]: Something went wrong with Diamond."
		exit 1 
	fi
fi
############################################################################################################################################################

if [[ $diamond -eq 1 ]]; then
### Relative abundances by mapping 
	printf '\n%s\n\n' "[INFO]: Counting abundances for Krona."

	bwa-mem2 index "$fasta"
	bwa-mem2 mem "$fasta" "$read1_path" "$read2_path" -t "$threads" | samtools view -Su - | samtools sort - -o "$sample".R.sort.bam
	bwa-mem2 mem "$fasta" "$unpaired_path" -t "$threads" | samtools view -Su - | samtools sort - -o "$sample".un.sort.bam
	samtools merge -f "$sample".bam "$sample".R.sort.bam "$sample".un.sort.bam
	rm "$sample".R.sort.bam
	rm "$sample".un.sort.bam
	samtools index "$sample".bam
	samtools idxstats "$sample".bam | cut -f1,3 > "$sample".magnitudes

### Krona visualization
	printf '\n%s\n\n' "[INFO]: Making Krona chart."
	ktImportBLAST -o "$sample".html "$sample".m8,"$sample" "$sample".m8:"$sample".magnitudes,"$sample".magn

	#ktClassifyBLAST "$outdir"/DIAMOND/"$sample".m8 -o KRONA/"$sample".tab
	#awk 'NR==FNR { a[$1]=$2; next} $1 in a {print $0,"\t"a[$1]}' "$outdir"/SCAFFOLDS/"$sample".magnitudes "$outdir"/KRONA/"$sample".tab > "$outdir"/KRONA/"$sample".magnitudes.tab
fi

if [[ $? -eq 0 ]]; then
	printf '\n%s\n' "[INFO]: viper-classify finished successfully! "
else
	>&2 printf '\n%s\n' "[ERROR]: viper-classify finished abnormally."
fi