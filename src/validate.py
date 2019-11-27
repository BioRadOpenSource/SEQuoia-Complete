"""
Short script to validate a fastq file. Will throw assertion error if any condition is not met
"""
import sys
import gzip

def validate_records(fastq_path):
    num_reads = 0
    with gzip.open(fastq_path, 'rb') as fh:
        for readname in fh:
            num_reads += 1
            seq = next(fh)
            sep = next(fh)
            qual = next(fh)

            assert len(seq) == len(qual)
            assert sep == b"+\n"
            assert readname.startswith(b'@')
    return num_reads

def main():
    read1 = sys.argv[1]
    num_read1 = validate_records(read1)
    assert num_read1 <= 500000000

    if len(sys.argv) == 3:
        read2 = sys.argv[2]
        num_read2 = validate_records(read1)
        assert num_read1 == num_read2

    print("Reads are valid", file=sys.stderr)


if __name__ == "__main__":
    main()
