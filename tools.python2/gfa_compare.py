#!/usr/bin/env python
"""
Compare two gfa files, reporting differences but ignoring some trivial ones
---------------------------------------------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

import sys,re

def usage(s=None):
	message = """
gfa_compare [--sort] gfa_file1 gfa_file2
"""

	if (s == None): sys.exit (message)
	else:           sys.exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	sortEm = False
	if (sys.argv[1] == "--sort"):
		sortEm = True
		del sys.argv[1]

	if (len(sys.argv) < 3):
		usage("you must specify two gfa files")
	elif (len(sys.argv) > 3):
		usage("wrong number of arguments")

	gfa1Filename = sys.argv[1]
	gfa2Filename = sys.argv[2]

	# compare the files

	gfa1 = file(gfa1Filename,"rt")
	gfa2 = file(gfa2Filename,"rt")

	if (sortEm): (different,lineNum) = compare_unsorted_files(gfa1,gfa2)
	else:        (different,lineNum) = compare_sorted_files  (gfa1,gfa2)

	if (different):
		print >>sys.stderr,"FAILURE: %s and %s are different (line %d)" \
		                 % (gfa1Filename,gfa2Filename,lineNum)
		sys.exit(1)

	print >>sys.stderr,"SUCCESS: %s and %s are equivalent" \
					 % (gfa1Filename,gfa2Filename)


# compare files that we expect are in the same order

def compare_sorted_files(gfa1,gfa2):
	lineNum = 0
	while (True):
		lineNum += 1
		line1 = gfa1.readline()
		line2 = gfa2.readline()
		if (line1 == "") and (line2 == ""):
			return (False,lineNum)
		line1 = line1.rstrip()
		line2 = line2.rstrip()

		stanza  = line1.split()[0]
		stanza2 = line2.split()[0]
		if (stanza2 != stanza):
			return (True,lineNum)

		if (stanza == "d"):
			continue	# ignore command line differences

		elif (stanza == "h"):
			line1 = " ".join(header_strip(line1))
			line2 = " ".join(header_strip(line2))

		if (line1 != line2):
			# print >>sys.stderr,"%s\n%s" % (line1,line2)
			return (True,lineNum)


# compare files that we suspect might not be in the same order

def compare_unsorted_files(gfa1,gfa2):

	lines1 = read_lines(gfa1)
	lines1.sort()

	lines2 = read_lines(gfa2)
	lines2.sort()

	compareNum = 0
	while (True):
		if (compareNum < len(lines1)): (line1,lineNum1) = lines1[compareNum]
		else:                          (line1,lineNum1) = ("",None)
		if (compareNum < len(lines2)): (line2,lineNum2) = lines2[compareNum]
		else:                          (line2,lineNum2) = ("",None)
		compareNum += 1
		if (line1 == "") and (line2 == ""):
			return (False,compareNum)

		stanza  = line1.split()[0]
		stanza2 = line2.split()[0]
		if (stanza2 != stanza):
			return (True,compareNum)

		if (stanza == "d"):
			continue	# ignore command line differences

		elif (stanza == "h"):
			line1 = " ".join(header_strip(line1))
			line2 = " ".join(header_strip(line2))

		if (line1 != line2):
			# print >>sys.stderr,"%s\n%s" % (line1,line2)
			return (True,compareNum)


# read all lines from a file, a return a list of (line,lineNumber)

def read_lines(gfa):
	lines = []
	lineNum = 0
	while (True):
		lineNum += 1
		line = gfa.readline()
		if (line == ""):
			break
		lines += [(line.rstrip(),lineNum)]
	return lines


headerRe = re.compile("^(?P<stanza>.+) +\"(?P<name1>.+)\" +\"(?P<name2>.+)\"$")

def header_strip(s):
	m = headerRe.match(s)
	if (m == None):
		return s

	stanza = m.group("stanza")
	name1  = m.group("name1").strip()
	name2  = m.group("name2").strip()

	if (name1.startswith(">")): name1 = name1[1:].strip()
	if (name2.startswith(">")): name2 = name2[1:].strip()

	return [stanza,name1,name2]


if __name__ == "__main__": main()
