#Script takes as input several graphics and stores the summary to RNASeq_internal_report.html

import sys, getopt, os, shutil, datetime
import argparse
import re


### parse options
parser = argparse.ArgumentParser(description="Creates a html report, storing all necessary files into html_summary. To visualize the report, keep the structure of the folder.")
parser.add_argument("-p", "--project", dest="project", help="name of the project")
parser.add_argument("-o", "--output", dest="output_dir", help="output directory")
parser.add_argument("-c", "--customer", dest="customer", help="group name of the customer", nargs="*")
parser.add_argument("-e", "--email", dest="email", help="email of the person who analyzed the data")
parser.add_argument("-n", "--names", dest="names", help="vector of all sample names", nargs="*")
parser.add_argument("-f", "--fastqc", dest="fastqc", help="vector of all sample names for the fastqc quality reports", nargs="*")
parser.add_argument("-u", "--unpaired", action="store_true", dest="unpaired", default=False, help="set if reads are unpaired (default: paired)")
parser.add_argument("-q", "--quality", dest="quality", help="per base sequence quality graph(s) (fastqc): <Sample1 R1> <Sample1 R2> <Sample2 R1> <Sample2 R2> [optional]", nargs="*")
parser.add_argument("-a", "--adapter", dest="adapter", help="adapter content graph(s) (fastqc): <Sample1 R1> <Sample1 R2> <Sample2 R1> <Sample2 R2> [optional]", nargs="*" )
parser.add_argument("-g", "--gccont", dest="gccontent", help="GC-content graph(s) (fastqc) [optional]", nargs="*")
parser.add_argument("-s", "--stat", dest="statistics", help="log file with alignment statistics (hisat2) [optional]", nargs="*")
parser.add_argument("-d", "--dist", dest="distribution", help="read distribution statistics text files (Rseqc, readDistributionStats.txt) [optional]", nargs="*")
parser.add_argument("-b", "--geneb", dest="genebody", help="gene body coverage graph (Rseqc) [optional]", nargs="*")
parser.add_argument("-i", "--insert", dest="insertsize", help="insert size histogram graph (rseqc) [optional]", nargs="*")
parser.add_argument("-t", "--sat", dest="saturation", help="saturation report graph(s) (Rseqc) [optional]", nargs="*")
parser.add_argument("--counts", dest="counts", help="table of counts with header (featureCounts) [optional]")
parser.add_argument("-y", "--type", dest="genetype", help="counts per gene biotype graph(s) [optional]", nargs="*")
parser.add_argument("-z", "--pca", dest="pca", help="PCA plot (DESeq2) [optional]")
parser.add_argument("-x", "--additional", dest="additional", help="Additional graph(s) you would like to append [optional]", nargs="*")


options = parser.parse_args()



### ------------------------ Functions -------------------------- ###


### create a labelled graph: the labelling is adjusted depending if the reads are unpaired or paired
def labelled_figure(options_argument, options_names):
	HTML_ELEMENT = ""
	if options.unpaired == False:
		l = len(options_names)
		for i in range(0,l):
			item1 = options_argument[(i+1)*2-2]
			item2 = options_argument[(i+1)*2-1]
			f1 = options_names[i] + "_R1_" + item1.split("/")[-1] #copy only filename without relative path
			f2 = options_names[i] + "_R2_" + item2.split("/")[-1]
			fig_name1 = "Sample " + options_names[i] + " R1"
			fig_name2 = "Sample " + options_names[i] + " R2"
			if item1 == "empty":
				s1 = (  '<figure class="block"><figure><figcaption>' + fig_name1 + '</figcaption><p>--- EMPTY PLOT ---</p></figure>\n\n')
			if item2 == "empty":
				s2 = (  '<figure class="block"><figure><figcaption>' + fig_name2 + '</figcaption><p>--- EMPTY PLOT ---</p></figure>\n\n')
			else:
				shutil.copy(item1, data_dir + f1)
				shutil.copy(item2, data_dir + f2)
				s1 = (  '<figure class="block"><figure><figcaption>' + fig_name1 + '</figcaption><img src="' + html_data + f1 + '" alt="' + fig_name1 + '"></figure>\n\n')
				s2 = (  '<figure><figcaption>' + fig_name2 + '</figcaption><img src="' + html_data + f2 + '" alt="' + fig_name2 + '"></figure></figure>\n\n')
			HTML_ELEMENT += s1
			HTML_ELEMENT += s2
	else:
		l = len(options_names)
		for i in range(0,l):
			item = options_argument[i]
			f = item.split("/")[-1] + options_names[i] #copy only filename without relative path
			fig_name = "Sample " + options_names[i]
			if item == "empty":
				s = (  '<figure><figcaption>' + fig_name + '</figcaption>><p>--- EMPTY PLOT ---</p></figure>\n\n')
			else:
				shutil.copy(item, data_dir + f)
				s = (  '<figure><figcaption>' + fig_name + '</figcaption><img src="' + html_data + f + '" alt="' + fig_name + '"></figure>\n\n')
			HTML_ELEMENT += s
	return HTML_ELEMENT

### create a simple figure
def simple_figure(options_argument):
	graph = options_argument
	f = graph.split("/")[-1]
	shutil.copy(graph, data_dir + f)
	if(is_img(options_argument)):
		s = '<img src="' + html_data + f + '" alt="' + f + '">\n\n'
	else:
		s = '<embed src="' + html_data + f + '" alt="' + f + '">\n\n'
	return s

### check format of input: return True if jpg or png, else False
def is_img(options_argument):
	last = options_argument[-3:]
	if(last == "jpg" or last == "png"):
		return True
	elif(last == "pdf"):
		return False
	else:
		print("\n\nWarning:\nImage format not jpg, png or pdf.\n\n")
		return False

### create simple labelled figure
def simp_fig(options_argument,label):
	graph = options_argument
	f = label + graph.split("/")[-1]
	shutil.copy(graph, data_dir + f)
	if(is_img(options_argument)):
		fig = '<figure><figcaption>' + label + '</figcaption><img src="' + html_data + f + '" alt="' + label + '"></figure>\n\n'
	else:
		fig = '<figure><figcaption>' + label + '</figcaption><embed src="' + html_data + f + '" alt="' + label + '"></figure>\n\n'
	return fig



def simple_labelled_figure(options_arguments):
	i = 0
	figure = ""
	for graph in options_arguments:
		label = "Sample " + options.names[i]
		figure += simp_fig(graph, label)
		i += 1
	return figure



#################################################################################
#################################################################################




### store output and graphs in folder output
output_dir = options.output_dir
html_data = "html_data_report/"
data_dir = os.path.join(output_dir, html_data)

if not os.path.exists(data_dir):
	os.makedirs(data_dir)



### create the content of the html file
now = datetime.datetime.now()
HTML_INDEX = '''
	<html>\n
		<head>\n
			<meta charset="utf-8">
			<title>RNA-Seq report ''' + options.project + '''</title>\n
			<style>
				body {font-family:Tahoma, Geneva, sans-serif;}
				h1 {font-size:16px;}
				h2 {font-size:14px}
				h3, p, figcaption, table {font-size:12px;}
				img {max-width: 500px;}
				figure {margin: 2px; display:inline-block;}
				.block {min-width: 1010px; display:block;}
				table {display:inline-block; text-align:left; border-collapse:collapse; margin-right:5px;margin-left:0px;}
				td {border: 1px solid black; border-collapse: collapse; padding:5px;}
				th {border: 1px solid black; border-collapse: collapse; padding:5px;font-weight:bold;}
				caption {text-align:left; padding-top:10px; padding-bottom: 5px;}
				embed {width:500px;min-height:500px;}
			</style>\n
		</head>\n
		<body>\n
			<h1>RNA-seq report: Project ''' + options.project + '''</h1>\n

		<p>
		Group: ''' + " ".join(options.customer) + '''</br>
		Report created: '''+ now.strftime("%d %B %Y") + '''</br>
		Analysis by: ''' + options.email +'''</p>\n\n'''


CLOSING_TAG = '</body></html>'




#################################################################################
#################################################################################




### 1. Read quality report
READ_QUAL = '\n<h2>1. Read quality report</h2>\n'
SEQ_QUAL = "\n\n"
ADP_CONT = "\n\n"
GC_CONT = "\n\n"

### 1.1 Per base sequence quality (fastqc)
if options.quality:
	SEQ_QUAL = '''<h3>1.1 Per base sequence quality (fastqc)</h3>\n
			<p>Base quality (Phred scores) along the length of the reads in each FastQ file. The box plots are drawn as follows: red line=median; yellow box=range between upper and lower quartiles; whiskers=range between 10% and 90% quantiles. The blue line shows the mean quality.</p>\n\n'''
	SEQ_QUAL += labelled_figure(options.quality, options.fastqc)


### 1.2 Adapter content (fastqc)
if options.adapter:
	ADP_CONT = '''<h3>1.2 Adapter content (fastqc)</h3>\n
			<p>Cumulative proportion of observed adapter sequences along the length of the reads in each FastQ file. This should be (close to) zero.</p>\n\n'''
	ADP_CONT += labelled_figure(options.adapter, options.fastqc)


### 1.3 GC content (fastqc)
if options.gccontent:
	GC_CONT = '''<h3>1.3 GC content</h3>\n
			<p>Distribution of GC content of the reads for each sample. We would expect a roughly normal distribution centred on the average GC content of the genome which varies between species.</p>\n\n'''
	GC_CONT += labelled_figure(options.gccontent, options.fastqc)


READ_QUAL += SEQ_QUAL + ADP_CONT + GC_CONT


### 2. Alignment report
ALIGN_REP = '\n<h2>2. Alignment report</h2>\n'
ALIGN_STAT = "\n\n"
READ_DIST = "\n\n"
GENE_BOD = "\n\n"
INS_SIZE = "\n\n"


### 2.1. Alignment statistics
### => Percent of reads mapped, total reads, unmapped reads
### For PE reads, ii)-iv) is for CONCORDANTLY mapped pairs only

if options.statistics and options.counts:
# It is enough to check if statistics and counts as a whole exits. There cannot be missing elements.
# Pipeline would crash earlier if there were.
# The order of the samples must be the same in names, counts and statistics as they are expanded in the same way
	ALIGN_STAT = '''<h3>2.1. Alignment statistics</h3>\n
		<p>Summary statistics for each sample:<br />
		i) Total = Total number of reads (for SE data) / read pairs (for PE data)<br />
		ii) Unique alignments = Number (percent) of reads / read pairs that can be assigned to exactly 1 position in the genome <sup>*)</sup><br />
		iii) Non-unique alignments = Number (percent) of reads / read pairs that match multiple positions in the reference genome <sup>*)</sup><br />
		iv) Unmapped = Number (percent) of reads / read pairs that cannot be mapped to the reference genome <sup>*)</sup><br />
		v) Total alignment rate (%)<br />
		vi) Aligned to genes = Number (percent) of reads overlapping with an annotated gene<br /><br />
		<sup>*)</sup> Note that, for PE data, these stats consider the number of pairs that map CONCORDANTLY 0, 1 or >1 times</p>\n\n'''

	def percent(number, total):
		perc = float(number)*100/total
		perc = str(round(perc,2))
		return  " (" + perc + "%) "

	row1 = "<tr><th></th>"
	row2 = "<tr><td>Total</td>"
	row3 = "<tr><td>Unique alignments</td>"
	row4 = "<tr><td>Non-unique alignments</td>"
	row5 = "<tr><td>Unmapped</td>"
	row6 = "<tr><td>Total alignment rate</td>"
	row7 = "<tr><td>Aligned to genes</td>"

	for i in range(0, len(options.names)):
		total_reads = 0
		total_rate = 'NA'
		unmapped_original = 0
		unique_original = 0
		multi_original = 0
		stats = open(options.statistics[i],"r")
		for line in stats:

			# values that are always the same:
			if re.search('reads; of these:', line):
				total_reads_original = int(line.split("r")[0].strip())
				total_reads = "{:,}".format(total_reads_original).replace(",", "'")
			elif re.search('overall alignment rate$', line):
				total_rate = line.split("o")[0]

			# values that depend on SE or PE
			else:
				if options.unpaired == False:
					if re.search('aligned concordantly 0 times$', line):
						unmapped_original = int(line.split("(")[0].strip())
						unmapped = "{:,}".format(unmapped_original).replace(",", "'")
					elif re.search('aligned concordantly exactly 1 time$', line):
						unique_original = int(line.split("(")[0].strip())
						unique = "{:,}".format(unique_original).replace(",", "'")
					elif re.search('aligned concordantly >1 times$', line):
						multi_original = int(line.split("(")[0].strip())
						multi = "{:,}".format(multi_original).replace(",", "'")
				else:
					if re.search('aligned 0 times$', line):
						unmapped_original = int(line.split("(")[0].strip())
						unmapped = "{:,}".format(unmapped_original).replace(",", "'")
					elif re.search('aligned exactly 1 time$', line):
						unique_original = int(line.split("(")[0].strip())
						unique = "{:,}".format(unique_original).replace(",", "'")
					elif re.search('aligned >1 times$', line):
						multi_original = int(line.split("(")[0].strip())
						multi = "{:,}".format(multi_original).replace(",", "'")
		stats.close()

### Obtain total number of reads aligned to genes from count table --> row 7
		countfile = open(options.counts,"r")
		genic_reads_original = 0
		for line in countfile:
			genic_reads_original += int(line.split("\t")[i+1])  # first column contains gene id
		genic_reads = "{:,}".format(genic_reads_original).replace(",", "'")
		countfile.close()

		row1 += "<th>Sample " + options.names[i] + "</th>"
		row2 += "<td>" + str(total_reads) + "</td>"
		row3 += "<td>" + str(unique) + percent(unique_original,total_reads_original) + "</td>"
		row4 += "<td>" + str(multi) + percent(multi_original,total_reads_original) + "</td>"
		row5 += "<td>" + str(unmapped) + percent(unmapped_original,total_reads_original) + "</td>"
		row6 += "<td>" + str(total_rate) + "</td>"
		row7 += "<td>" + str(genic_reads) + percent(genic_reads_original,total_reads_original) + "</td>"

	row1 += "</tr>"
	row2 += "</tr>"
	row3 += "</tr>"
	row4 += "</tr>"
	row5 += "</tr>"
	row6 += "</tr>"
	row7 += "</tr>"

	table = "\n\n<table>\n" + row1 + row2 + row3 + row4 + row5 + row6 + row7 + "\n</table>\n\n"
	ALIGN_STAT += table


### 2.2 Read distribution
if options.distribution:
	READ_DIST = '''<h3>2.2 Read distribution</h3>\n
			<p>Overview of how the reads are distributed among different genomic features (TSS=Transcription start site; TES=Transcription end site). Total_bases=Total length of each category in the genome; Tag_count=Number of reads assigned to each category. A read may be split and assigned to multiple categories; Tags/kb=Number of reads per kilobase of genomic sequence.</p>\n\n'''
	j = 0
	for file in options.distribution:
		table = '\n\n<table>\n<caption>Sample ' + options.names[j] + '</caption>'
		if file == "empty":
			table += '\n<tr><th></th><th></th><th></th><th></th></tr><tr><td></td><td></td><td></td><td></td></tr></table>\n\n'
		else:
			read_dist = open(file, "r")
			read_dist = read_dist.readlines()

			row = read_dist[4].split()
			th = '<tr><th>'+row[0]+'</th>'+'<th>'+row[1]+'</th>'+'<th>'+row[2]+'</th>'+'<th>'+row[3]+'</th></tr>'
			table += th

			for i in range(5,len(read_dist)-1):
				row = read_dist[i].split()
				tr = '<tr><td>'+row[0]+'</td>'+'<td>'+format(int(row[1]),',d')+'</td>'+'<td>'+format(int(row[2]),',d')+'</td>'+'<td>'+row[3]+'</td></tr>'
				table = table + tr
		table += '\n</table>\n\n'
		READ_DIST += table
		j += 1


### 2.3 Gene body coverage
if options.genebody:
	GENE_BOD = '''<h3>2.3 Gene body coverage</h3>\n
			<p>Distribution of reads along the length of the genes (5' end on left, 3' on right)</p>\n\n'''
	GENE_BOD += simple_labelled_figure(options.genebody)


### 2.4 Insert size histogram
if options.insertsize:
	INS_SIZE = '''<h3>2.4 Insert size histogram</h3>\n
			<p>Histogram of inner distance for each sample (available for PE data only). Inner distance is the length of the RNA fragment minus the length of the two reads.</p>\n\n'''
	INS_SIZE += simple_labelled_figure(options.insertsize)


ALIGN_REP += ALIGN_STAT + READ_DIST + GENE_BOD + INS_SIZE




### 3. Counts report_REP = '<h2>3. Counts report</h2>'
COUNTS_REP = '<h2>3. Counts report</h2>\n'
PCA_PLOT = "\n\n"
SAT_REP = "\n\n"
COUNTS_PER_GENE = "\n\n"

### 3.1 PCA
if options.pca:
	PCA_PLOT = '''<h3>3.1 PCA</h3>\n
			<p>Plot of the first two axes from a principal component analysis based on the 500 genes with the most variable expression across all samples. Different colours indicate the different experimental groups.</p>\n\n'''
	PCA_PLOT += simple_figure(options.pca)



### 3.2 Saturation report
if options.saturation:
	SAT_REP = '''<h3>3.2 Saturation report</h3>\n
			<p>The number of splice junctions detected using different subsets of the data from 5% to 100% of all reads. Red=known junction based on the provided genome annotation; green=novel junctions; blue=all. At sequencing depths sufficient to perform alternative splicing analyses, at least the red line should reach a plateau where adding more data does not much increase the number of detected junctions.</p>\n\n'''
	SAT_REP += simple_labelled_figure(options.saturation)



### 3.3 Counts per gene type
if options.genetype:
	COUNTS_PER_GENE = '''<h3>3.3 Counts per gene type</h3>\n
				<p>Boxplots of the distribution of counts in each sample across all genes within a particular category. Only genes with count > 0 are considered. Note that this step needs to download info from biomart which sometimes does not work. In this case, an empty plot will be displayed.</p>\n\n'''
	COUNTS_PER_GENE += simple_labelled_figure(options.genetype)



COUNTS_REP += PCA_PLOT + SAT_REP + COUNTS_PER_GENE

### write html
HTML_INDEX += READ_QUAL + ALIGN_REP + COUNTS_REP + CLOSING_TAG
outputfile = open(os.path.join(output_dir, "RNASeq_internal_report.html"), "w+")


outputfile.write(HTML_INDEX)
