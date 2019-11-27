import fileinput
from collections import Counter


def main():
    start, end = 0, 8
    for line in fileinput.input():
        bad = False
        vals = line.rstrip("\n").split("\t")
        barcode = vals[5][start:end]
        barcode_len = len(barcode)
        missing = 8 - barcode_len
        # If the length - number of N's is less than 7, throw it out
        # basically allows an edit distance of 1
        if missing + Counter(barcode).get("N", 0) > 1:
            bad = True
            barcode = "N" * (end - start)

        vals[0] = "@{}_{}".format(barcode, vals[0][1:])
        vals[4] = "@{}_{}".format(barcode, vals[4][1:])
        if not bad:
            print("\t".join(vals))
        else:
            print("MM{}".format("\t".join(vals)))


if __name__ == "__main__":
    main()
