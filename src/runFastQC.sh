    
#!/usr/bin/env bash
set -e

#BEGIN ARGUMENT PARSING
function usage
{
	echo "	";
    echo "usage: docker run -v /path/to/sample/:/data/ bioraddbg/fastqc -i input_dir -o output_dir";
    echo "	";
    echo "	-i | --input_directory		: Input directory containing fastq files for QC";
    echo "	-o | --output_directory		: Output directory to write files";
    echo "	-h | --help			: This message";
    echo "	";
}

function parse_args
{
	# positional args
	args=()

	# named args
	while [ "$1" != "" ]; do
		case "$1" in
			-i | --input_directory )		input_dir="$2";		shift;;
			-o | --output_directory )		output_dir="$2";		shift;;
			-h | --help )					usage;					exit;; # quit and show usage
			* )								args+=("$1")			# if no match, add it to the positional args
		esac
	shift # move to next kv pair
	done

	# restore positional args
	set -- "${args[@]}"

	#Input data directory selected?
	if [[ -z "$input_dir" ]]; then
		echo "No input directory selected."
		usage
		exit 1;
	fi

	#Input data directory exists?
	if [[ ! -d "$input_dir" ]]; then
		echo "Input directory does not exist -- provide a directory with fastq files for QC"
		exit 1;
	fi

	#Input data has fastq files?
	if [[ -z "$(find "$input_dir" -maxdepth 1 | grep fastq)" ]]; then
		echo "Input directory has no FASTQ files.";
		exit 1;
	fi
}

function run
{
  parse_args "$@"

  echo "you passed in..."
  echo "Input Directory: $input_dir"
  echo "Output Directory: $output_dir"
}

run "$@";

#BEGIN FUNCTION

set -x

# input_dir="$1"
# output_dir="$2"

mkdir -p "$output_dir"
fastqs=$(find "$input_dir" -maxdepth 1 -type f -name '*fastq*')
for fastq in $fastqs
do
	sem --bg -j $(nproc) fastqc $fastq -O "$output_dir"
done
sem --wait

