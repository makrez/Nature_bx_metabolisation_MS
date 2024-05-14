# The order of the samples is the same in names, counts and statistics as they are expanded in the same way


import argparse, re

parser = argparse.ArgumentParser(description="Creates a table with mapping statistics for customer report.")
parser.add_argument("-n", "--names", dest="names", help="vector of all sample names", nargs="*")
parser.add_argument("-s", "--stat", dest="statistics", help="log file with alignment statistics (hisat2)", nargs="*")
parser.add_argument("--counts", dest="counts", help="table of counts with header (featureCounts) [optional]")
parser.add_argument("-u", "--unpaired", action="store_true", dest="unpaired", default=False, help="set if reads are unpaired (default: paired)")
parser.add_argument("-o", "--output", dest="output", help="output file")
options = parser.parse_args()



header = ["Sample", "Nbr reads", "% mapped", "Nbr multi-mapped", "Nbr genic"]
header = "\t".join(header)
outputfile = open(options.output, "w")
outputfile.write(header + "\n")


for i in range(0, len(options.names)):

	stats = open(options.statistics[i], "r")

	temp = [options.names[i], "NA", "NA", "NA", "NA"]

	for line in stats:
# values that are always the same:
		if re.search('reads; of these:', line):
			total_reads_original = int(line.split("r")[0].strip())
			total_reads_temp = "{:,}".format(total_reads_original).replace(",", "'")
		elif re.search('overall alignment rate$', line):
			total_rate_temp = line.split("o")[0]

# values that depend on SE or PE
		else:
			if options.unpaired == False:
				if re.search('aligned concordantly >1 times$', line):
					multi_original = int(line.split("(")[0].strip())
					multi_temp = "{:,}".format(multi_original).replace(",", "'")
			else:
				if re.search('aligned >1 times$', line):
					multi_original = int(line.split("(")[0].strip())
					multi_temp = "{:,}".format(multi_original).replace(",", "'")
	stats.close()

	temp[1] = total_reads_temp
	temp[2] = total_rate_temp
	temp[3] = multi_temp

### Obtain total number of reads aligned to genes from count table --> row 7
	countfile = open(options.counts,"r")
	genic_reads_original = 0
	for line in countfile:
		genic_reads_original += int(line.split("\t")[i+1])  # first column contains gene id
	genic_reads_temp = "{:,}".format(genic_reads_original).replace(",", "'")
	countfile.close()

	temp[4] = genic_reads_temp

	temp = "\t".join(str(element) for element in temp)
	outputfile.write(temp + "\n")

outputfile.close()
