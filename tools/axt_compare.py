#!/usr/bin/env python
"""
Compare two axt files, reporting differences but ignoring some trivial ones
---------------------------------------------------------------------------

:Author: Bob Harris (rsharris@bx.psu.edu)
"""

import sys,re

def usage(s=None):
	message = """
axt_compare [--sort] axt_file1 axt_file2
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
		usage("you must specify two axt files")
	elif (len(sys.argv) > 3):
		usage("wrong number of arguments")

	axt1Filename = sys.argv[1]
	axt2Filename = sys.argv[2]

	# compare the files

	axt1 = file(axt1Filename,"rt")
	axt2 = file(axt2Filename,"rt")
	different = compare_files(axt1,axt2,sortEm=sortEm)
	axt1.close()
	axt2.close()

	if (type(different) == tuple):
		(message,which) = different
		if   (which == "axt1"): message = "in %s, %s" % (axt1Filename,message)
		elif (which == "axt2"): message = "in %s, %s" % (axt2Filename,message)
		different = True
	elif (type(different) == str):
		message = "%s and %s are different, %s" \
		        % (axt1Filename,axt2Filename, different)
		different = True
	elif (different):
		message = "%s and %s are different" \
		        % (axt1Filename,axt2Filename)

	if (different):
		print >>sys.stderr,"FAILURE: %s" % message
		sys.exit(1)

	print >>sys.stderr,"SUCCESS: %s and %s are equivalent" \
					 % (axt1Filename,axt2Filename)


# compare files

def compare_files(axt1,axt2,sortEm=False):

	blocks1 = read_axt_blocks(axt1)
	if (type(blocks1) == str): return (blocks1,"axt1")

	blocks2 = read_axt_blocks(axt2)
	if (type(blocks2) == str): return (blocks2,"axt2")

	if (len(blocks1) != len(blocks2)):
		return ("different number of blocks",None)

	if (sortEm):
		blocks1 = sort_blocks(blocks1)
		blocks2 = sort_blocks(blocks2)

	for ((block1,lineNum1),(block2,lineNum2)) in zip(blocks1,blocks2):
		if (type(block1[0]) == str): block1 = convert_block(block1)
		if (type(block2[0]) == str): block2 = convert_block(block2)
		info1 = block1[0]
		if (len(info1) != 9): return ("bad axt block at line %d" % lineNum1,"axt1")
		info2 = block2[0]
		if (len(info2) != 9): return ("bad axt block at line %d" % lineNum2,"axt2")

		if (sortEm):
			block1[0][0] = 0 # (clear axt record number)
			block2[0][0] = 0 # (clear axt record number)

		if (block1 == block2): continue
		return "block at line %d vs block at line %d" % (lineNum1,lineNum2)

	return False


# sort the blocks in a list, by everything but the axt record number

def sort_blocks(blocks):
	newBlocks = []
	for (block,lineNum) in blocks:
		(info,text1,text2) = newBlock = convert_block(block)
		key = info[1:]
		newBlocks += [(key,newBlock,lineNum)]
	newBlocks.sort()
	return [(block,lineNum) for (key,block,lineNum) in newBlocks]


# convert a block from textual representation to something more amenable to
# comparisons and sorting

def convert_block(block):
	(info,text1,text2) = block

	info = info.split()
	return (info,text1,text2)


# read all blocks from a file, a return a list of (block,lineNumber)

def read_axt_blocks(f):
	blocks = []
	block = None
	lineNum = 0
	for line in f:
		lineNum += 1
		line = line.strip()
		if (line == "") or (line.startswith("#")):
			if (block != None):
				if (len(block) != 3):
					return ("bad axt block at line %d" % blockLineNum)
				blocks += [(block,blockLineNum)]
			block = None
			continue
		if (block == None):
			block = []
			blockLineNum = lineNum
		block += [line]
	if (block != None):
		if (len(block) != 3):
			return ("bad axt block at line %d" % blockLineNum)
		blocks += [(block,blockLineNum)]
	return blocks


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
