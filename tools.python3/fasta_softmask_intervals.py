#!/usr/bin/env python3
"""
Given a list of intervals, mask those bases in the fasta sequence(s).
"""

from sys import argv,stdin,exit


def usage(s=None):
	message = """fasta_softmask_intervals [options] < fasta_file > fasta_file
  Apply masking intervals to create a soft-masked fasta file.

  options:
    <intervals_file>          file containing a list of intervals to be masked,
                              in the form <chrom> <start> <end>;  --origin
                              determines whether these are origin one or zero
    --complement              mask the intervals *not* in the intervals file
                              (not to be confused with reverse-complementation)
    --chrom=<sequence_names>  (cumulative) copy (and mask) only the specified
                              sequence(s) <sequence_names> is a comma-separated
                              list
                              (default is to copy and mask all sequences)
    --origin=one              intervals are origin-one, closed
                              (default is origin-zero, half-open)
    --wrap=<line_length>      split each sequence into multiple lines if needed
                              (default is to write sequence on a single line)
    --mask=<character>        mask with a particular character (usually X or N)
                              (default is to mask with lowercase)"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse args

	chromsOfInterest = None
	doComplement     = False
	origin           = "zero"
	wrapLength       = 100
	maskChar         = None
	intervalsFile    = None

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg in ["--help","-h","--h","-help"]):
			usage()
		elif (arg.startswith("--chrom=")) or (arg.startswith("--chroms=")):
			if (chromsOfInterest == None):
				chromsOfInterest = []
			chromsOfInterest += argVal.split(",")
		elif (arg == "--complement"):
			doComplement = True
		elif (arg.startswith("--origin=")):
			origin = argVal
			if (origin == "0"): origin = "zero"
			if (origin == "1"): origin = "one"
			if (origin not in ["zero","one"]):
				usage("unknown argument: %s=%s" % (arg,val))
		elif (arg.startswith("--wrap=")):
			wrapLength = int(argVal)
		elif (arg.startswith("--mask=")):
			maskChar = argVal
			if (len(maskChar) != 1): usage("--mask requires a single character")
		elif (arg.startswith("--")):
			usage("can't understand %s" % arg)
		elif (intervalsFile == None):
			intervalsFile = arg
		else:
			usage("can't understand %s" % arg)

	if (intervalsFile == None):
		usage("you have to tell me the intervals you're interested in")

	# read the intervals

	f = open(intervalsFile,"rt")

	chromToIntervals = {}

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == "") or (line.startswith("#")): continue

		fields = line.split()
		assert (len(fields) >= 3), \
		      "not enough fields (line %s): %s" % (lineNumber,line)

		try:
			chrom  = fields[0]
			start = int(fields[1])
			end   = int(fields[2])
			if (origin == "one"): start -= 1
			if (start < 0):    raise ValueError
			if (start >= end): raise ValueError
		except ValueError:
			assert (False), \
			      "bad line (line %s): %s" % (lineNumber,line)

		if (chromsOfInterest != None) and (chrom not in chromsOfInterest):
			continue

		if (chrom not in chromToIntervals): chromToIntervals[chrom] = []
		chromToIntervals[chrom] += [(start,end)]

	f.close()

	for chrom in chromToIntervals:
		chromToIntervals[chrom] = merge_and_sort(chromToIntervals[chrom])

	# process the sequences

	chromSeen = {}

	for (chrom,seq) in fasta_sequences(stdin):
		if (chromsOfInterest != None) and (chrom not in chromsOfInterest):
			continue

		assert (chrom not in chromSeen), \
			"more than one sequence is named %s" % chrom
		chromSeen[chrom] = True

		seq = seq.upper()
		if (chrom not in chromToIntervals): chromToIntervals[chrom] = []

		intervals = chromToIntervals[chrom]
		if (doComplement): intervals = complement_intervals(intervals,len(seq))

		newSeq = []

		prevEnd = 0
		for (start,end) in intervals:
			if (prevEnd < start):  newSeq += [seq[prevEnd:start]]
			if (maskChar == None): newSeq += [seq[start:end].lower()]
			else:                  newSeq += [maskChar*(end-start)]
			prevEnd = end
		if (prevEnd < len(seq)):   newSeq += [seq[prevEnd:]]

		print (">%s" % chrom)
		newSeq = "".join(newSeq)
		assert (len(newSeq) == len(seq)), "internal error"

		for i in range(0,len(newSeq),wrapLength):
			print ("".join(newSeq[i:i+wrapLength]))

	# make sure all sequences were given

	missing = [chrom for chrom in chromToIntervals if (chrom not in chromSeen)]
	assert (missing == []), "missing fasta sequence %s" % (", ".join(missing))


# fasta_sequences--
#	Read the fasta sequences from a file

def fasta_sequences(f):
	seqName = None
	seqNucs = None

	for line in f:
		line = line.strip()

		if (line.startswith(">")):
			if (seqName != None):
				yield (seqName,"".join(seqNucs))
			seqName = line[1:].strip().split()[0]
			seqNucs = []
		elif (seqName == None):
			assert (False), "first sequence has no header"
		else:
			seqNucs += [line]

	if (seqName != None):
		yield (seqName,"".join(seqNucs))


# merge_and_sort--
#	Merge a set of intervals (union of sets) and sort them by increasing
#	position

def merge_and_sort(intervals):
	intervals.sort()

	start = None
	for (s,e) in intervals:
		if (start == None):
			(start,end) = (s,e)
		elif (s > end):
			yield (start,end)
			(start,end) = (s,e)
			continue
		elif (e > end):
			end = e

	if (start != None):
		yield (start,end)


# complement_intervals--
#	Produce the set-wise complement of a set of intervals. We assume the
#	intervals are sorted by increasing position, and that there are no overlaps.

def complement_intervals(intervals,seqLen):
	prevEnd = 0
	for (start,end) in intervals:
		if (start != prevEnd):
			yield (prevEnd,start)
		prevEnd = end
	if (seqLen != prevEnd):
		yield (prevEnd,seqLen)


if __name__ == "__main__": main()
