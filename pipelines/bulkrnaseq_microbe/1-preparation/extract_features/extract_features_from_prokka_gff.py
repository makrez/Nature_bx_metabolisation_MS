'''
This script extracts the following info from a prokka generated gff file:
gene id, gene symbol, product name
the resulting list can be added to RNAseq DEG table of
microbial RNAseq DE analysis pipeline results in the case a
custom reference is used and therefor Biomart can not be used
How to run the script:
python extract_genesymbol_from_gff.py --gff foo.gff

You will need python3.7 to run this script!

Simone Oberhaensli, 2021-08-31
'''


import argparse
import re
import csv




def main():
    # input is the gene annotation in gff format
    parser = argparse.ArgumentParser(description='description')
    parser.add_argument("--gff", help="a gff file")
    args = parser.parse_args()

    # read gff
    with open(args.gff, "r") as gff:
        lines = gff.readlines()
    #''' loop through gff line by line, extract lines which contain information on
    #gene ID, symbol and product and extract these info from the line. Store info in a list
    #and print it to a file in csv format
    #'''
    genes = []
    # define these variables according to your gff
    ptr0 = 'CDS' # tRNAs, ribosomal RNAs are missed
    ptr1 = ';gene='
    ptr2 = ';product='
    ptr3 = 'ID='
    for line in lines:
        if re.search(ptr0, line):
            #id = re.search(rf'{ptr3}(\S+);', line)
            id = re.search(rf'{ptr3}(\S+);Parent', line)
            symbol = re.search(rf'{ptr1}(.+);i', line) # is followed by inference
            product = re.search(rf'{ptr2}(.+)', line) # is followed by nothing
            #apparently some genes don't have a symbol
            if symbol is not None:
                lst = [id.group(1),symbol.group(1),product.group(1)]
                genes.append(lst)
            else:
                lst = [id.group(1), 'NA', product.group(1)]
                genes.append(lst)

    # create output list
    with open("features.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerows(genes)



if __name__ == "__main__":
    main()

