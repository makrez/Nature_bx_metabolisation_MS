import os, re, glob, sys, string
import pandas as pd

# make lists of all samples and all file names and a dictionary of all files that belong to a given sample
def getInfo(fastqfolder, extension, delim, start, end, mate1):
	FilesPerSample={}
	SampleList=[]
	FileBasename=[]
	uniqueSamples=[]

	for file in os.listdir(fastqfolder):
		if file.endswith(extension):

			temp=file.split(delim)[(start-1):end]
			temp=delim.join(temp) #ok also when len(temp) is 1
			temp=temp+delim # makes sure that full sample ID is matched
			SampleList.append(temp)

			# Make a list of all file basenames.
			if file.endswith(delim+mate1+"."+extension):
				temp2=re.sub(mate1+"."+extension, "", file)
				FileBasename.append(temp2)

	uniqueSamples=sorted(list(set(SampleList)), key=str.lower)
	# Make a dictionary that has sample IDs as keys and lists of all file basenames that belong to that sample as values
	for x in uniqueSamples:
		temp3=[]
		for y in FileBasename:
			if y.startswith(x):  # sample name will always be at the beginning of the file name
				temp3.append(y)
		FilesPerSample[x]=temp3

	FileBasename=sorted(list(set(FileBasename)), key=str.lower)

	return [FileBasename, uniqueSamples, FilesPerSample]

### Function to set up vector of all contrasts
# Takes all factor levels from experimental design file and returns a vector of all
# possible pairwise combinations
# If there is a level called control, it will be put as the second element in the contrast to ensure that
# log foldchanges will be output as group X / control. Case is ignored
def Contrasts(expDesign):
	import itertools
	inputfile=open(expDesign, "r")
	levels = []
	for line in inputfile:
		temp1=line.strip('\n').strip('\r')
		if temp1=='':
			continue
		else:
			levels.append(str(temp1.split('\t')[1]))
	inputfile.close()
	levels = sorted(list(set(levels)))
	# If there is a level called Control, make sure it is put last
	temp = [x.lower() for x in levels]
	try:
		c_index = temp.index('control') # temp is only used to find the position
		levels.append(levels.pop(c_index))
	except ValueError:
		levels = levels

	contrasts = []
	for comb in itertools.combinations(levels, 2):
		contrasts.append(".".join(comb))
	return(contrasts)

def readContrasts(contrastfile):
	contrasts = []
	with open(contrastfile, "r") as infile:
		for line in infile:
			combs = line.strip().split('\t')
			if(combs == ['']):
				continue
			if(len(combs) != 2):
				sys.exit('Check your contrasts file. Each row should contain exactly two tab separated columns.')
			comb = ".".join(combs)
			contrasts.append(comb)
	return(contrasts)

# 1. Check that the sample names are the same in expDesign and in fastq files:
# create a vector containing all sample IDs that are listed in the expDesign file
def getExpDesignIDs(expDesign, delim):
	expDesignIDs = []
	design = open(expDesign, 'r')

	for line in design:
		row = line.strip('\n').strip('\r')
		row = row.split('\t')
		row = row[0]
		# check that ID ends with delim. Add if not
		if not row == '':
			if not row.endswith(delim):
				row = row + delim
			expDesignIDs.append(row)

	expDesignIDs = sorted(expDesignIDs, key=str.lower)
	return(expDesignIDs)

# Create a vector containing all sample IDs from the file names

def getfileIDs(fastqFolder, extension, delim, start, end):
	fileIDs = []
	fastqFiles = glob.glob(fastqFolder + "/*." + extension)

	fastqFiles = [x.replace(str(fastqFolder), "") for x in fastqFiles]

	for file in fastqFiles:
		file = file.split(delim)
		file = file[start-1:end]
		file = delim.join(file)
		fileIDs.append(file + delim)
	fileIDs = sorted(set(fileIDs), key=str.lower)
	return(fileIDs)

def checkFileExists(Files):
	for x in Files:
		if not os.path.isfile(x):
			sys.exit("File {} does not exist".format(x))

def checkRpath(Rpath):
	if not os.path.isdir(Rpath):
		sys.exit("R path is not correct")

def checkHisat2Indexes(reference):
	if (len(glob.glob(reference + ".*ht2"))!=8 and len(glob.glob(reference + ".*ht2l"))!=8):
		sys.exit("Check that all Hisat2 index files exist! There should be 8 files ending in .[1-8].ht2 or .[1-8].ht2l")

# 2. Check that the sampleName is correct. It cannot contain any special characters other than delim AND
# it cannot start with a number

def checkSpecialCharacters(expDesignIDs):
	allowed = string.ascii_uppercase + string.ascii_lowercase + "".join([str(e) for e in range(0,10)]) + "_"

	for id in expDesignIDs:
		if all(c in allowed for c in id) and not id[0].isdigit():
			continue
		else:
			sys.exit("At least one of your sample IDs contains a special character or starts with a number. Sample IDs must start with a letter (lower or uppercase) and contain only letters, numbers and underscores '_'")

# Convert memory from GB to MB. Drmaa recommendation is to use mem_mb.
# In the web form, the format will be checked to make sure it is as expected (class memory)
def convert_mem(mem_gb):
	mem_mb = int(mem_gb.replace('G', ''))*1000
	return mem_mb

# For each kit, output the correct strandedness setting
def infer_strandedness(kit):
	strandinfo = {
	"Truseq": 2,
	"Corall": 1,
	"Takara": 2,
	"Quantseq": 0
	}
	try:
		return(strandinfo[kit])
	except Exception as e:
		sys.exit("This is not a valid kit name.")



# Check number of columns in exp design file and summarize groups and replicates
def checkDesign(design_file, model):
	df = pd.read_table(design_file, header = None)

	# Check number of columns:
	# exp design file has to contain 2 columns for uni-factorial analysis
	# 3 columns for two-factorial
	columns_expected = {
		"unifactorial": 2,
		"twofactorial": 3
	}
	if not (columns_expected[model] == len(df.columns)):
		sys.exit("The number of columns in your experimental design file does not match the model you specified. For a uni-factorial design there must be two columns (SampleID, ExperimentalGroups). For a two-factorial design, there must be three columns (SampleID, ExperimentalGroups, Confounder)." + "\n\n" + "You selected a " +  model + " model. Your experimental design file contains " + str(len(df.columns)) + " columns.")

	# Print tables with the number of samples per group
	print("Your experimental design file specifies these experimental groups and number of replicates per group:")
	if(model == "unifactorial"):
		print(df.iloc[:, 1].value_counts())
	else:
		df.columns = ["SampleID", "Groups", "Confounder"]
		print(df[["Groups","Confounder"]].value_counts())



