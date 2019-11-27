""" Test Data
Geneid	Length	ExpectedTPM	ExpectedRPKM	TotalReads
A	2000	333333.333333	142857.142857	10
B	4000	333333.333333	142857.142857	20
C	1000	333333.333333	142857.142857	5
D	10000	0	0	0
"""

import pandas as pd
import sys


def main():
    """Add a RPKM and TPM Column to the input file.

    RPKM - Reads Per Kilobase Million (read_counts / (total_reads/1_000_000)) / gene_length_kb
    TPM  - Transcripts Per Kilobase Million (read_counts / gene_length_kb) / (sum(read_counts / gene_length_kb) / 1_000_000)
    """
    scale = 1000000
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    count_data = pd.read_csv(input_file, sep="\t", comment="#")
    count_data.rename(columns={count_data.columns[-1]: "TotalReads"}, inplace=True)

    # RPKM
    total_reads = count_data["TotalReads"].sum()
    total_reads_per_mil = total_reads / scale
    count_data["RPKM"] = count_data.apply(
        lambda x: (x["TotalReads"] / total_reads_per_mil) / (x["Length"] / 1000), axis=1
    )

    # TPM
    rpk = count_data.apply(lambda x: x["TotalReads"] / (x["Length"] / 1000), axis=1)
    sum_rpk_per_mil = rpk.sum() / scale
    count_data["TPM"] = rpk.apply(lambda x: x / sum_rpk_per_mil)
    count_data.to_csv(path_or_buf=output_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
