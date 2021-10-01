#!/usr/bin/env python
"""
Given pairwise alignments from lastz's general format, perform "chaining"
(essentially syntenic filtering) on them.
"""

from sys    import argv,stdin,stdout,stderr,exit
from math   import ceil


def usage(s=None):
	message = """

WARNING: THIS MODULE HAS NOT BEEN TESTED. USE AT YOUR OWN PERIL.

usage: cat alignments | chain_alignments [options]
  --format=<list>      provide comma-separated list of the names of the
                       columns, in order; these must include the field names
                       that are listed in detail below
  --format=auto        read column names from the first line of the input,
                       which must begin with a "#"
  --chain=<diag,anti>  penalties for diagonal and anti-diagonal
                       (default penalties are 0)
  --match=<reward>     score for matched bases; this informs us of the scale
                       of the scores given in the input; in lastz this was
                       taken from the scoring matrix's A-to-A match value
                       (default value is 91)

The input is alignments in lastz's general format. Required columns are name1,
start1, end1, name2, strand2, start2, end2, and score. Note that strand1 is
present it must be "+". The length of the alignment in each species must be the
same. Columns can appear in any order, and other columns may be included.
Alignments need not be sorted in any particular order. A typical example is
shown below.

  #name1 start1 end1  name2  start2 end2  strand2 score
  APPLE  10436  10727 ORANGE 1      292   +       19786
  APPLE  2112   2404  ORANGE 1      293   +       19230
  APPLE  1      242   ORANGE 291    532   +       15055
  APPLE  11775  11999 ORANGE 538    762   +       13045
  APPLE  11533  11767 ORANGE 523    757   +       13506
   ...

Output is in the same format as input. Alignments for the same group (same 
name pair and strand) will appear together (i.e. on consecutive lines). No
guarantee is made on the order of alignments within a group.

WARNING: THIS MODULE HAS NOT BEEN TESTED. USE AT YOUR OWN PERIL."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


requiredColumns = ["name1","start1","end1",
                   "name2","strand2","start2","end2",
                   "score"]
optionalColumns = ["strand1"]
columnAliases   = {"s"  : "strand2",
                   "s2" : "strand2"}


def main():
	global columnNames,headerLine
	global debug

	chainScale  = 100.0  # (this is hardwired and unchangeable in lastz)

	# parse the command line

	columnNames = None
	diagPenalty = 0.0
	antiPenalty = 0.0
	aaMatch     = 91.0   # score for A-A match
	debug       = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg == "--format=auto") or (arg == "--format=automatic"):
			columnNames = "automatic"
		elif (arg.startswith("--format=general:")):
			argVal = argVal.split(":",1)[1]
			columnNames = argVal.split(",")
		elif (arg.startswith("--format=")):
			columnNames = argVal.split(",")
		elif (arg == "--chain"):
			diagPenalty = antiPenalty = 0.0
		elif (arg.startswith("G=")):
			diagPenalty = float(argVal)
		elif (arg.startswith("R=")):
			antiPenalty = float(argVal)
		elif (arg.startswith("--chain=")):
			(diagPenalty,antiPenalty) = argVal.split(",",1)
			diagPenalty = float(diagPenalty)
			antiPenalty = float(antiPenalty)
		elif (arg.startswith("--match=")):
			aaMatch = float(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (columnNames == None):
		usage("you must tell me the input column names")
	elif (columnNames == "automatic"):
		columnNames = None
	else:
		columnNames = column_names(columnNames)

	print >>stderr, "WARNING: THIS MODULE HAS NOT BEEN TESTED. USE AT YOUR OWN PERIL."

	# collect alignments by species pair and strand

	pairStrandToAlignments = {}
	pairsSeen = set()
	pairs = []

	for a in read_alignments(stdin):
		pair = (a.name1,a.name2)
		if (pair not in pairsSeen):
			pairsSeen.add(pair)
			pairs += [pair]

		pairStrand = (pair,a.strand2)
		if (pairStrand not in pairStrandToAlignments):
			pairStrandToAlignments[pairStrand] =  [a]
		else:
			pairStrandToAlignments[pairStrand] += [a]

	# output chained alignments; the order of species pairs will be the same
	# as observed in the input, except that we output each pair's + strand
	# before its - strand

	headerPrinted = False

	for pair in pairs:
		for strand2 in ["+","-"]:
			pairStrand = (pair,strand2)
			if (pairStrand not in pairStrandToAlignments):
				continue

			(best,chain) = reduce_to_chain \
			                 (pairStrandToAlignments[pairStrand],
			                  diagPenalty,antiPenalty,chainScale,aaMatch,
			                  chain_connect_penalty)

			for a in chain:
				if (not headerPrinted) and (headerLine != None):
					print headerLine
					headerPrinted = True

				print a.line


# read_alignments--

class Alignment: pass

def read_alignments(f):
	global headerLine
	global columnNames
	global name1Column,strand1Column,start1Column,end1Column
	global name2Column,strand2Column,start2Column,end2Column
	global scoreColumn

	headerLine    = None
	columnsNeeded = None

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line.startswith("#")):
			headerLine = line
			if (columnNames == None):
				fields = line.split()
				fields[0] = fields[0][1:]
				columnNames = column_names(fields)
				idColumn = None
			continue

		assert (columnNames != None), \
		       "input column names are not provided within the file"

		if (columnsNeeded == None):
			columnsNeeded = 1 + max([columnNames[name] for name in columnNames])

			(name1Column,start1Column,end1Column,
			 name2Column,strand2Column,start2Column,end2Column,
			 scoreColumn) = requiredColumns
			(strand1Column,) = optionalColumns

			name1Column   = columnNames[name1Column  ]
			strand1Column = columnNames[strand1Column] if (strand1Column in columnNames) else None
			start1Column  = columnNames[start1Column ]
			end1Column    = columnNames[end1Column   ]
			name2Column   = columnNames[name2Column  ]
			strand2Column = columnNames[strand2Column]
			start2Column  = columnNames[start2Column ]
			end2Column    = columnNames[end2Column   ]
			scoreColumn   = columnNames[scoreColumn ]

		fields = line.split()
		assert (len(fields) >= columnsNeeded), \
		       "not enough columns at line %d (%d, expected %d)" \
		     % (lineNumber,len(fields),columnsNeeded)

		a = Alignment()
		a.lineNumber = lineNumber
		a.line       = line
		a.name1      = fields[name1Column  ]
		a.strand1    = fields[strand1Column] if (strand1Column != None) else None
		a.pos1       = fields[start1Column ]
		end1         = fields[end1Column   ]
		a.name2      = fields[name2Column  ]
		a.strand2    = fields[strand2Column]
		a.pos2       = fields[start2Column ]
		end2         = fields[end2Column   ]
		a.score      = fields[scoreColumn ]

		try:
			a.pos1 = int(a.pos1)
			end1   = int(end1) + 1
			if (a.pos1 >= end1): raise ValueError
			a.length = end1-a.pos1
		except ValueError:
			assert (False), \
			       "bad alignment (at line %d), first species start/end\n%s" \
			     % (lineNumber,line)

		assert (a.strand1 in [None,"+","-"]), \
		       "bad alignment (at line %d), first species strand\n%s" \
		     % (lineNumber,line)
		assert (a.strand1 != "-"), \
		       "bad alignment (at line %d), first species strand must be \"+\"\n%s" \
		     % (a.lineNumber,a.line)

		try:
			a.pos2 = int(a.pos2)
			end2   = int(end2) + 1
			if (a.pos2 >= end2): raise ValueError
		except ValueError:
			assert (False), "bad alignment (at line %d), second species start/end\n%s" \
			              % (lineNumber,line)

		assert (a.strand2 in ["+","-"]), \
		       "bad alignment (at line %d), second species strand\n%s" \
		     % (lineNumber,line)

		assert (end2-a.pos2 == a.length), \
		       "bad alignment (at line %d), unequal lengths\n%s" \
		     % (lineNumber,line)

		try:
			a.score = float(a.score)
			if (a.score <= 0): raise ValueError
		except ValueError:
			assert (False), \
			       "bad alignment (at line %d), bad score\n%s" \
			     % (lineNumber,line)

		yield a


# column_names--

def column_names(names):
	columnNames = {}
	for (ix,name) in enumerate(names):
		actualName = name
		if (name in columnAliases): name = columnAliases[name]
		if (name not in requiredColumns + optionalColumns): continue
		if (name in columnNames):
			usage("\"%s\" (or alias) appears more than once in --format" % actualName)
		columnNames[name] = ix
	for name in requiredColumns:
		if (name not in columnNames):
			usage("--format lacks required name \"%s\"" % name)
	return columnNames


#----------
#
# chain--
#	Find the highest scoring chain in a set of gap-free alignments. Each
#	segment in the chain will begin strictly before the start of the next
#	segment. This is (expected to be) the most parsimonious subset of the
#	gap-free alignments, assuming there actual orthology contains no
#	inversions.
#
# The algorithm finds, for each segment, the highest scoring chain that ends
# with that segment. Segments are scanned in an order (by increasing start in
# sequence 1) that guarantees all possible predecessor chains have been found
# and scored before that segment is considered. Upon completion, the chain is
# recovered by backtracking from its end segment.
#
# A chain's score is the sum of its segment scores minus the sum of penalties
# for the gaps between segments. The caller must provide a function to compute
# those penalties. See note (1) of reduce_to_chain() for more details.
#
# To facilitate the search for valid predecessors, a K-d tree is used. See the
# header of build_kd_tree() for more details on the tree implementation.
#
# References:
#
#	[1] Multidimensional Binary Search Trees Used for Associative Searching.
#	    Jon Louis Bentley, Commun. ACM 18(9): 509-517 (1975).
#
#----------

# reduce_to_chain--
#	Find the highest scoring chain, in which each segment in the chain begins
#	strictly before the start of the next segment.
#
#	This implementation is intended to be a faithful copy of reduce_to_chain()
#	from lastz's chain.c. However, there are a few differences, primarily in
#	the calling interface.
#
# A chain is a series of segments, where each segment in the chain (other than
# the last), begins strictly before the start of the next. A chain's score is
# scale times the sum of segment scores minus the sum of penalties for the gaps
# between segments:
#		connect (segment_i, segment_(i+1), scale)
# the last sum is taken over all segments in the chain except the last).
#
# Arguments:
#	alignments:				The alignment segments on which to operate.
#	diagPenalty:			Chaining penalty;  see notes (1) and (3).
#	antiPenalty:			Chaining penalty;  see notes (1) and (3).
#	scale:					Scaling constant;  see note (2).
#	aaMatch:				Alignment score of an A-to-A match.
#	connectionPenalizer:	Chain connection penalty function;  see note above,
#							.. and description of arguments in chain.h
#
# Returns:
#	A pair, (best,alignments):
#		best:		The score of the best chain, unscaled;  zero if there's
#					.. some problem.
#		alignments:	The alignments belonging to the chain.
#
# Notes:
#	(1)	The parameters diagPenalty and antiPenalty permit us to deduce useful
#		inequalities about chain scores. Namely, let segment_i and segment_j
#		be segments on diagonals diag_i and diag_j, and set
#			diff = diag_j - diag_i
#		Then diagPenalty and antiPenalty are required to satisfy:
#			if diff >= 0, then connect(segment_i,segment_j) >= diff*diagPenalty
#		and
#			if diff < 0, then connect(segment_i,segment_j) >= -diff*antiPenalty
#
#	(2)	In effect, scale permits integer arithmetic to be used with very small
#		gap penalties, since the computed chain also maximizes the sum of the
#		segment scores minus the sum of
#			connect(segment_i, segment_(i+1), scale)/scale.
#
#	(3) diagPenalty and antiPenalty are considered to have already been scaled.
#		We only apply scale to the segment substitution scores.

class ChainingState: pass

class BestPred: pass

def reduce_to_chain(alignments,
                    diagPenalty,antiPenalty,scale,aaMatch,
                    connectionPenalizer):
	global state

	n = len(alignments)
	if (n == 0): return (0.0,[])

	state = ChainingState()
	state.diagPenalty = diagPenalty
	state.antiPenalty = antiPenalty
	state.scale       = scale
	state.aaMatch     = aaMatch
	state.connect     = connectionPenalizer

	# sort segments by pos1, so that the predecessor search loop is guaranteed
	# to score all possible predecessors of any segment before it considers
	# that segment

	segments = [(seg.pos1,seg) for seg in alignments]
	segments.sort()
	state.segments = [seg for (_,seg) in segments]

	state.perm       = [0] * n
	state.invPerm    = [0] * n
	state.chainScore = [0.0] * n

	# build the K-d tree; as part of this process the segments are permuted
	# (by use of the perm[] array), and we compute the inverse of that
	# permutation to aid later access to the segments

	state.perm = list(xrange(n))			# build the identity permutation,
											# .. mapping node numbers to
											# .. segments
	root = build_kd_tree(0,n-1,1)			# build the K-d tree (alters
											# .. state.perm)
	for i in xrange(n):						# compute the inverse permutation
		 state.invPerm[state.perm[i]] = i

	# for each segment, find the best chain ending at that segment; the array
	# chain[] provides the path from any segment back through its chain; any
	# segment i for which best_predecessor() finds no positive scoring
	# predecessor (such as those that have no predecesssor at all) will have
	# chain[i]==None and terminate backtracking

	bp = BestPred()
	chain = [0] * n

	best    = 0
	bestEnd = None
	for i in xrange(n):
		state.query = state.segments[i]
		state.x     = state.query.pos1
		state.y     = state.query.pos2
		state.diag  = state.x - state.y

		bp.num     = None
		bp.contrib = 0
		best_predecessor(bp,root,1,0)  # (modifies bp)
		queryContrib = state.query.score * state.scale
		state.chainScore[i] = queryContrib + bp.contrib

		if (state.chainScore[i] > best):
			best = state.chainScore[i]
			bestEnd = i
		chain[i] = bp.num
		propagate_max_score(root,state.chainScore[i],state.invPerm[i])

	# backtrack to collect chain segments

	keep = [False] * n

	i = bestEnd
	while (i != None):
		keep[i] = True
		i = chain[i]

	chainedSegments = [state.segments[i] for i in xrange(n) if (keep[i])]

	# scale back best score and return the chainedSegments

	return (best/scale,chainedSegments)


# chain_connect_penalty--
#	Compute penalty for connecting two segments in the chain.
#
#	This implementation is intended to be a faithful copy of
#	chain_connect_penalty() from lastz's lastz.c.

def chain_connect_penalty(seg1,seg2):
	# nota bene: numSubs is the number of substitutions needed to get from end
	#            .. of segment 1 to beginning of segment 2

	assert ((seg2.pos1 > seg1.pos1) and (seg2.pos2 > seg1.pos2)), \
	       "HSPs improperly ordered for chaining\n%s\n%s" \
	     % (seg1.line,seg2.line)

	xEnd  = seg1.pos1 + seg1.length - 1
	yEnd  = seg1.pos2 + seg1.length - 1

	diag1 = seg1.pos1 - seg1.pos2
	diag2 = seg2.pos1 - seg2.pos2

	diagDiff = diag2 - diag1
	if (diagDiff >= 0):
		# segment 1's diagonal is above segment 2's
		numSubs = seg2.pos2- yEnd - 1
	else:
		# segment 1's diagonal is below segment 2's
		numSubs  = seg2.pos1 - xEnd - 1
		diagDiff = -diagDiff

	penalty = diagDiff * state.diagPenalty
	if (numSubs >= 0):
		penalty +=   numSubs  * state.antiPenalty
	else:
		penalty += (-numSubs) * state.scale * state.aaMatch

	# nota bene: in lastz, the penalty was clipped to he maximum possible
	#            score (0x7FFFFFFF for integer scoring) to prevent overflow;
	#            we don't perform such clipping in this implementation

	return penalty


# build_kd_tree--
#	Build segments into a K-d tree.
#
# Standard K-d tree implimentation (for K=2), such as might be found in
# reference [1]. The points are partitioned into two sets, split by the a value
# along one axis. Each of those sets is in turn split again, along the other
# axis, and so on, until all sets are small enough. "Small enough" is defined
# by bucketSize. The two dimensional axes are y (sequence pos2) and diagonal
# (sequence pos1-pos2).
#
# Arguments:
#	lo,hi:	range of entries (of state.segments[], indexed by state.perm[]) to
#			.. build a tree of; these are inclusive (i.e. there are hi+1-lo
#			.. entries)
#	axis:	which dimension/axis to partition (at the top level)
#			  0 => diagonal (pos1 - pos2)
#			  1 => pos2
#
# Returns:
#	The root of the tree. state.perm[] is modified so that entries in
#	state.segments[state.perm[]] agree with the tree.

class KDNode: pass
#	isBucket				# True => this node is a bucket/leaf
#	loIx, hiIx				# if isBucket is True
#							# ..   index range of the segments in this leaf
#							# if isBucket is False
#							# ..   hiIx is index corresponding to cutVal
#	# the following fields are only valid if isBucket is False
#	cutVal					# signed value (along appropriate axis) which
#							# .. separates lower and upper children
#	maxChainScore			# the highest score for any chain ending at a
#							# .. segment in this subtree
#	loSon,hiSon				# pointers to child nodes

bucketSize = 3				# max number of entries we'll place in a bucket node


def build_kd_tree(lo,hi,axis):
	p = KDNode()
	p.maxChainScore = 0

	if (hi+1-lo <= bucketSize):		# the range is small enough to fit in one
		p.isBucket = True			# .. node
		p.loIx   = lo
		p.hiIx   = hi
		p.cutVal = None
		p.loSon  = None
		p.hiSon  = None
	else:							# the range is too big for one node, split
		p.isBucket = False			# .. it into two subtrees
		m = partition_segments(lo,hi,axis)
		p.loIx   = None
		p.hiIx   = m
		p.cutVal = projection(m,axis)
		p.loSon  = build_kd_tree(lo, m, 1-axis)
		p.hiSon  = build_kd_tree(m+1,hi,1-axis)

	assert (p != None), \
	       "(in build_kd_tree, p == None)"

	assert (valid_KDNode(p)), \
	       "(in build_kd_tree, p is not a valid KDNode)"

	return p


# partition_segments--
#	Partition a list of segments into two sets, split by a value along one axis.
#
# A 'pivot' value is desginated as the median (along the specified axis) of
# the first and last segment, and the one at the middle of the list. Segments
# below the pivot are moved to the first part of the list, and segments above
# it are moved to the last part, with the pivot between them.
#
# Arguments:
#	lo,hi:	range of entries (of state.segments[], indexed by state.perm[]) to
#			.. build a tree of; these are inclusive (i.e. there are hi+1-lo
#			.. entries)
#	axis:	dimension/axis to partition on
#			  0 => diagonal (pos1 - pos2)
#			  1 => pos2
#
# Returns:
#	The index of the pivot (m). state.perm[] is modified so that entries in
#	state.segments[state.perm[]] satsify lo..m-1 <= m <= m+1..hi.

def partition_segments(lo,hi,axis):
	assert (hi - lo >= 2), \
	       "partition: cannot happen (%d,%d)" \
	     % (lo, hi)

	while (True):
		# find the pivot and move it to the front; we use the median of the
		# lower, middle, and upper values as the pivot

		m = (lo+hi) // 2
		a = projection(lo,axis)
		b = projection(m, axis)
		c = projection(hi,axis)

		if ((a <= b) and (b <= c)) or ((c <= b) and (b <= a)):
			perm_swap(lo,m)
			pivot = b
		elif ((a <= c) and (c <= b)) or ((b <= c) and (c <= a)):
			perm_swap(lo,hi)
			pivot = c
		else:
			pivot = a

		# move smaller entries to front, larger to back

		i = lo
		j = hi+1
		while (i < j):
			# search forward for a large entry
			i += 1
			while (i<=hi) and (projection(i,axis)<=pivot):
				i += 1

			# search backward for a small entry
			j -= 1
			while (j>=lo) and (projection(j,axis)>pivot):
				j -= 1

			perm_swap(i,j)

		perm_swap(i,j)		# undo the last swap
		perm_swap(lo,j)		# move the pivot value to the proper location

		# warning: we must avoid returning j==hi (because build_kd_tree() would
		# recurse forever); if j<hi, we had at least one value larger than the
		# pivot; when j==hi the pivot was the max value; if the range had only
		# two values we are assured that these two are sorted

		if   (j < hi):     return j
		elif (hi-lo == 2): return hi-1

		# otherwise, we need to partition them again, leaving out hi; looping
		# back is equivalent to tail recursion, with this call:
		#    return partition_segments (lo, hi-1, axis)

		hi -= 1


# best_predecessor--
#	Find a segment's best predecessor chain.
#
# The best predecessor of a segment is the chain, among all those starting
# strictly before the segment in both x and y, that scores the highest when
# connected to this segment.
#
# This routine searches a subtree for such predecessor chains, pruning
# subtrees which contain no predecessor segments, or in which no segment can
# exceed the given lower bound.
#
# Arguments:
#	bp:			The best predecessor found so far. This will be modified in
#				.. place.
#	subtree:	The K-d (sub)tree to search.
#	axis:		Dimension/axis to partition on
#				  0 => diagonal (pos1 - pos2)
#				  1 => pos2
#	lowerBound:	Lower bound of chain score that must be achieved.
#
# Returns, by modifying bp:
#	The ending segment index (bestpred.num) and score (bestpred.contrib) of the
#	best predecessor chain.

def best_predecessor(bp,subtree,axis,lowerBound):
	assert (subtree != None), \
	       "(in best_predecessor, NULL subtree)"

	if (bp.contrib >= subtree.maxChainScore - lowerBound):
		return

	assert (valid_KDNode(subtree)), \
	       "(in best_predecessor, invalid subtree)"

	# if we're at a leaf, search over all segments in the leaf

	if (subtree.isBucket):
		for i in xrange(subtree.loIx,subtree.hiIx+1):
			j = state.perm[i]     # j is the segment we want to add to the chain
			s = state.segments[j] # s is the candidate to be a predecessor
			if (s.pos1 >= state.x) or (s.pos2 >= state.y):
				continue
			predScore = state.chainScore[j] - state.connect(s,state.query)
			if (predScore > bp.contrib):
				bp.contrib = predScore
				bp.num = j

	# if we're at a node cut by y, search over both subtrees, pruning the high
	# subtree if all its segments have y greater than our query

	elif (axis == 1):
		if (state.y >= subtree.cutVal):
			best_predecessor(bp,subtree.hiSon,lowerBound,1-axis)
		best_predecessor(bp,subtree.loSon,lowerBound,1-axis)

	# if we're at a node cut by the diagonal, search over search both subtrees,
	# adjusting the lower bound accordingly
	# nota bene: diff>0 => query diagonal is below cut

	else: # if (axis == 0)
		diff = state.diag - subtree.cutVal
		if (diff >= 0):	# query diagonal is southeast of (or same as) cut
			best_predecessor(bp,subtree.hiSon,1-axis,lowerBound)
			best_predecessor(bp,subtree.loSon,1-axis,diff*state.diagPenalty)
		else:			# query diagonal is northwest of cut
			best_predecessor(bp,subtree.loSon,1-axis,lowerBound)
			best_predecessor(bp,subtree.hiSon,1-axis,-diff*state.antiPenalty)


# propagate_max_score--
#	Propagate the best score for any chain ending at a particular segment to
#	all the (sub)trees that contain that segment.
#
# Arguments:
#	subtree:	The K-d (sub)tree to operate on.
#	s:			The score to propagate.
#	ix:			The index of the segment that has that score.

def propagate_max_score(subtree,s,ix):
	while (subtree != None):
		if (s > subtree.maxChainScore):
			subtree.maxChainScore = s
		if (ix <= subtree.hiIx): subtree = subtree.loSon
		else:                    subtree = subtree.hiSon


# projection-- figure out spatial position of segment i along the current axis

def projection(i,axis):
	if (axis == 0):
		return state.segments[state.perm[i]].pos1 - state.segments[state.perm[i]].pos2
	else:
		return state.segments[state.perm[i]].pos2


def perm_swap(i,j):
	(state.perm[i],state.perm[j]) = (state.perm[j],state.perm[i])


def valid_KDNode(p):
	return (p.isBucket) or ((p.loSon != None) and (p.hiSon != None))


if __name__ == "__main__": main()
