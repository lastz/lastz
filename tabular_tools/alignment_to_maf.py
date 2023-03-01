#!/usr/bin/env python3
"""
Convert pairwise alignments from lastz's general format to MAF.

References
  [1] MAF format (specification)
        http://genome.ucsc.edu/FAQ/FAQformat#format5
"""

from sys                  import argv,stdin,stderr,exit
from math                 import ceil
from alignment_table      import AlignmentTable
from alignment_table_aids import construct_alignment_text,score_alignment
from cigar                import split_cigar

try:                from two_bit_file import TwoBitFile
except ImportError: TwoBitFile = None

try:                from fasta_file import FastaFile
except ImportError: FastaFile = None


programName    = "alignment_to_maf"
programVersion = "0.1.0"


def usage(s=None):
	message = """
usage: cat <alignment_file> | %s [options]
  --sequence[s]=<file>    (required) specify the reference sequence(s); <file>
						  can be .2bit, .fasta, or .fasta.gz; for .fasta, if an
						  associated .fai file exists, it is utilized
  --alias:<alias>=<name>  the input can use <alias> as an alias for column
                          name <name>
  --head=<number>         limit the number of alignment records
  --progress=<number>     periodically report how many alignment records we've
                          processed
  --version               report this program's version number

Read pairwise alignments from lastz's general format and, using the cigar
strings, create a dotplot table suitable for plotting in R.

Typical input is shown below but other columns may be included, and columns can
appear in any order. Note that cigar or cigarx (either one) is required. If
score is included, it is used; otherwise the score is re-calculated (using
lastz defaults). If text1 and/or text2 are both included, they are used;
otherwise they are reconstructed using the cigar string. 

  #name1    zstart1  end1     name2            strand2 zstart2+  end2+     cigarx
  mule.chr1 1238615  1650901  donkey.contig12  +       1183966   1594899   ...
  mule.chr1 1952507  2309884  donkey.contig12  +       1889643   2247263   ...
  mule.chr1 3530276  3888853  donkey.contig12  +       88799552  89157262  ...
  mule.chr1 62122793 62459747 donkey.contig163 +       154012987 154349815 ...
  mule.chr1 4856000  5077740  donkey.contig163 +       90480442  90701655  ...
  mule.chr1 3926177  4143832  donkey.contig163 +       89190167  89407555  ...
   ...

Output is in MAF format, as described at
  http://genome.ucsc.edu/FAQ/FAQformat#format5""" \
% programName

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	fieldJoiner = "\t"

	# parse the command line

	referenceFilenames = []
	aliasToName        = {}
	headLimit          = None
	reportProgress     = None
	debug              = []
	refDebug           = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--sequence=")) or (arg.startswith("--sequences=")) or (arg.startswith("--ref=")):
			referenceFilename = argVal
			if  ((referenceFilename.endswith(".2bit")) \
			  or (referenceFilename.endswith(".fa.gz")) \
			  or (referenceFilename.endswith(".fasta.gz")) \
			  or (referenceFilename.endswith(".fsa_nt.gz")) \
			  or (referenceFilename.endswith(".fa")) \
			  or (referenceFilename.endswith(".fasta")) \
			  or (referenceFilename.endswith(".fsa_nt"))):
				pass
			else:
				usage("unrecognized file extension in %s" % arg)
			referenceFilenames += [referenceFilename]
		elif (arg.startswith("--alias:")):
			argVal = arg.split(":",1)[1]
			for argField in argVal.split(","):
				(alias,name) = argField.split("=",1)
				assert (alias not in aliasToName)
				aliasToName[alias] = name
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg in ["--version","--v","--V","-version","-v","-V"]):
			exit("%s, version %s" % (programName,programVersion))
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg in ["--refdebug","--2bitdebug"]):
			refDebug += ["debug"]
		elif (arg.startswith("--refdebug=")) or (arg.startswith("--2bitdebug=")):
			refDebug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (referenceFilenames == []):
		usage("you need to provide a reference genome file")

	if (len(referenceFilenames) > 1):
		usage("sorry, support for more than one reference genome file had not been implemented yet")

	# open the reference lookup (if we have one); this will map genomic
	# positions or segments to the corresponding nucleotide(s)

	refLookup = None

	if (referenceFilenames != []):
		referenceFilename = referenceFilenames[0]

		if (referenceFilename.endswith(".2bit")):
			if (TwoBitFile == None):
				assert (False), \
					   "two bit format support is not installed (needed for \"%s\")" % referenceFilename
			refLookup = TwoBitFile(referenceFilename,unmask=False,
		                           debug=refDebug)
		elif (referenceFilename.endswith(".gz")) or (referenceFilename.endswith(".gzip")):
			# (assume it is gzipped fasta)
			if (FastaFile == None):
				assert (False), \
					   "fasta format support is not installed (needed for \"%s\")" % referenceFilename
			refLookup = FastaFile(referenceFilename,unmask=False,
			                      ignoreFai=True,
			                      debug=refDebug)
		else:
			# (assume it is fasta)
			if (FastaFile == None):
				assert (False), \
					   "fasta format support is not installed (needed for \"%s\")" % referenceFilename
			refLookup = FastaFile(referenceFilename,unmask=False,
			                      ignoreFai=("ignorefai" in debug),
			                      debug=refDebug)

	# open the alignment table, as an iterator
	#
	# nota bene: we require the cigarx field but allow "cigar" as an alias for
	#            it; in this way we'll guarantee we have one or the other; and
	#            either is adequate for our purposes; note that this will be
	#            a.cigarx in an alignment object regardless of which the user
	#            provided

	requiredColumns    = ("name1","zstart1","end1",
	                      "name2","zstart2+","end2+",
	                      "strand",
	                      "cigarx")
	nonRequiredColumns = ("score","text1","text2")
	columnAliases      = {"strand2" : "strand",
	                      "s"       : "strand",
	                      "s2"      : "strand",
	                      "cigar"   : "cigarx"}
	for alias in aliasToName:
		columnAliases[alias] = aliasToName[alias]

	t = AlignmentTable.from_file(stdin,
	                             requiredColumns=requiredColumns,
	                             nonRequiredColumns=nonRequiredColumns,
	                             columnAliases=columnAliases)

	# process the alignments

	headerHasBeenWritten = False

	alignmentNumber = 0
	for a in t:
		alignmentNumber += 1

		if (headLimit != None) and (alignmentNumber > headLimit):
			print("limit of %s alignments reached" % ("{:,}".format(headLimit)),file=stderr)
			break
		if (reportProgress != None):
			if (alignmentNumber == 1) or (alignmentNumber % reportProgress == 0):
				print("processing alignment %s" % "{:,}".format(alignmentNumber),file=stderr)

		if (not headerHasBeenWritten):
			header = ["##maf","version=1"]
			header += ["scoring=lastz_defaults"]
			print(" ".join(header))
			print("# this file generated by %s version %s, converting from lastz tabular format" \
			     % (programName,programVersion))
			print("#")
			print("# gap_open_penalty   = 400")
			print("# gap_extend_penalty = 30")
			print("#        A    C    G    T")
			print("#   A   91 -114  -31 -123")
			print("#   C -114  100 -125  -31")
			print("#   G  -31 -125  100 -114")
			print("#   T -123  -31 -114   91")
			headerHasBeenWritten = True

		if ("linenumber" in debug):
			print("# line %d" % a.lineNumber)
		alignment_to_maf(a,refLookup)


# alignment_to_maf--
#	print a (pairwise) alignment in maf format

def alignment_to_maf(a,refLookup):
	seq1 = refLookup[a.name1]
	assert (seq1 != None), "no sequence has been provided for %s" % (a.name1)
	seq2 = refLookup[a.name2]
	assert (seq2 != None), "no sequence has been provided for %s" % (a.name2)
	(srcSize1,srcSize2) = (len(seq1),len(seq2))

	if (hasattr(a,"text1")) and (hasattr(a,"text2")):
		(text1,text2) = (a.text1,a.text2)
	else:
		cigarOps = split_cigar(a.cigarx)
		(text1,text2) = construct_alignment_text(a,cigarOps,refLookup,refLookup)

	if (hasattr(a,"score")):
		score = a.score
	else:
		score = score_alignment(text1,text2)

	start1 = a.start1
	size1  = a.end1 - a.start1

	start2 = a.start2 if (a.strand == "+") else (srcSize2 - a.end2)
	size2  = a.end2 - a.start2

	srcW     = max(len(a.name1)      ,len(a.name2))
	startW   = max(len(str(start1  )),len(str(start2  )))
	sizeW    = max(len(str(size1   )),len(str(size2   )))
	srcSizeW = max(len(str(srcSize1)),len(str(srcSize2)))

	print("a score=%s" % score)
	print("s %-*s %*s %*s %s %*s %s" \
	    % (srcW,a.name1,
	       startW,start1,
	       sizeW,size1,
	       "+",
	       srcSizeW,srcSize1,
	       text1))
	print("s %-*s %*s %*s %s %*s %s" \
	    % (srcW,a.name2,
	       startW,start2,
	       sizeW,size2,
	       a.strand,
	       srcSizeW,srcSize2,
	       text2))
	print("")


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return          int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))


if __name__ == "__main__": main()
