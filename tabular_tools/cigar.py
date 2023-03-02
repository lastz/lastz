#!/usr/bin/env python3
"""
Support for alignment cigar-strings.

References
  [1] Sequence Alignment/Map Format Specification (the CIGAR field in
      section 1.4)
        https://samtools.github.io/hts-specs/SAMv1.pdf
"""


# cigar_to_counts--
#	Report the number of match, mismatch, and indel columns in the alignment
#	represented by a cigar string.
#
#   Note that the cigar string should have "X" and "=" operators instead of
#	simply "M".

def cigar_to_counts(cigar):
	cigarInfo = split_cigar(cigar)
	if (cigarInfo == None): return None

	numMatch = numMismatch = numIndel = 0

	for (rpt,op) in cigarInfo.operations:
		if (op == "="):
			numMatch += rpt
		elif (op == "X"):
			numMismatch += rpt
		elif (op in ["I","D"]):
			numIndel += rpt
		elif (op == "M"):
			raise ValueError("unsupported operation \"%d%s\" in cigar \"%s\""
			               + "\nonly cigar strings with operators in {=,X,I,D} are supported"
			               % (rpt,op,cigar))
		else:
			raise ValueError("unsupported operation \"%d%s\" in cigar \"%s\""
			               % (rpt,op,cigar))

	return (numMatch,numMismatch,numIndel)


# cigar_to_windowed_counts--
#	Report (yield) the number of match, mismatch, and indel columns in the
#	alignment represented by a cigar string, yielding the count for each
#	overlapping window in the reference sequence.
#
#   Note that the cigar string should have "X" and "=" operators instead of
#	simply "M".
#
#   Per reference [1],
#     "I" is a insertion to the reference; it consumes the reference but not the
#         query
#     "D" is a deletion from the reference; it consumes the query but not the
#         reference

def cigar_to_windowed_counts(cigar,windowSize,windowStep=1):
	cigarInfo = split_cigar(cigar)
	if (cigarInfo == None): return None

	# convert the cigar operations to a string

	opStr = []
	for (rpt,op) in cigarInfo.operations:
		if (op in ["=","X","I","D"]):
			opStr += [op*rpt]
		elif (op == "M"):
			raise ValueError("unsupported operation \"%d%s\" in cigar \"%s\""
			               + "\nonly cigar strings with operators in {=,X,I,D} are supported"
			               % (rpt,op,cigar))
		else:
			raise ValueError("unsupported operation \"%d%s\" in cigar \"%s\""
			               % (rpt,op,cigar))

	opStr = "".join(opStr)

	# assign in-string positions to each reference location

	refLocToIx = []
	refLoc = 0
	for (ix,ch) in enumerate(opStr):
		if (ch in ["=","X","I"]):
			refLocToIx += [ix]   # nota bene: refLocToIx[refLoc] = ix
			refLoc += 1
		else: # if (ch == "D"):
			pass
	refLocToIx += [len(opStr)]   # nota bene: refLocToIx[refLoc] = len(opStr)

	refLen = len(refLocToIx) - 1

	# extract window-sized substrings and count the operations within; note
	# that we strip any prefix or suffix deletion since an aligner would not
	# include those in an alignment of the interval

	refLoc = 0
	while (refLoc+windowSize <= refLen):
		ixLft = refLocToIx[refLoc]
		ixRgt = refLocToIx[refLoc+windowSize]
		windowStr = opStr[ixLft:ixRgt].strip("D")

		numMatch    = windowStr.count("=")
		numMismatch = windowStr.count("X")
		numIndel    = windowStr.count("I") + windowStr.count("D")
		yield (refLoc,numMatch,numMismatch,numIndel)

		refLoc += windowStep


# Earlier implementation of cigar_to_windowed_counts. This might be expected to
# be more efficient than "B", at least in terms of memory used.
#
# This seems to be correct but with one exception. If a window has a deletion
# prefix, this are unintentionally counted as part of the window's indels; I
# wasn't able to see an easy way to correct that. It seems like any time I
# yield after consuming an operation, I would need to look ahead to see if the
# next operation is a deletion.

#def cigar_to_windowed_counts_A(cigar,windowSize):
#	cigarInfo = split_cigar(cigar)
#	if (cigarInfo == None): return None
#	numOps = len(cigarInfo.operations)
#
#	numMatch = numMismatch = numIndel = 0
#	lftIx,lftRpt = (-1,0)
#	rgtIx,rgtRpt = (-1,0)
#	currentLen = 0   # length of current window, in the reference/target sequence
#	while (True):
#		# if we've exhausted the right-op fetch the next one, and validate it
#
#		if (rgtRpt == 0):
#			rgtIx += 1
#			if (rgtIx == numOps): break
#			(rgtRpt,rgtOp) = cigarInfo.operations[rgtIx]
#			if (rgtOp == "M"):
#				raise ValueError("unsupported operation \"%d%s\" in cigar \"%s\""
#							   + "\nonly cigar strings with operators in {=,X,I,D} are supported"
#							   % (rgtRpt,rgtOp,cigar))
#			elif (rgtOp not in ["=","X","I","D"]):
#				raise ValueError("unsupported operation \"%d%s\" in cigar \"%s\""
#							   % (rgtRpt,rgtOp,cigar))
#
#		# if we haven't yet accumulated a full window, absorb the right-op
#		# (some or all of it) into the first window; if as a result the window
#		# is full, yield the counts; note that once we fill the window this
#		# if clause will never trigger again
#
#		if (currentLen < windowSize):
#			shortfall = windowSize - currentLen
#			if (rgtOp == "="):
#				rptAbsorbed = min(rgtRpt,shortfall)
#				numMatch   += rptAbsorbed
#				rgtRpt     -= rptAbsorbed
#				currentLen += rptAbsorbed
#			elif (rgtOp == "X"):
#				rptAbsorbed = min(rgtRpt,shortfall)
#				numMismatch += rptAbsorbed
#				rgtRpt      -= rptAbsorbed
#				currentLen  += rptAbsorbed
#			elif (rgtOp == "I"):
#				# consumes reference only
#				rptAbsorbed = min(rgtRpt,shortfall)
#				numIndel   += rptAbsorbed
#				rgtRpt     -= rptAbsorbed
#				currentLen += rptAbsorbed
#			else: # if (rgtOp == "D"):
#				# consumes query only
#				numIndel += rgtRpt
#				rgtRpt = 0
#
#			if (currentLen == windowSize):
#				yield (..,numMatch,numMismatch,numIndel)
#			continue
#
#		# otherwise, the window is full, and we'll absorb an equal amount of
#		# left- and right-ops (equal relative to the reference sequence); note
#		# that we don't track currentLen here, because an invariant is that it
#		# will always be equal to the windowSize
#
#		# if we've exhausted the left-op fetch the next one; note that we don't
#		# need to validate it because we know we did so earlier when it was a
#		# right-op
#
#		if (lftRpt == 0):
#			lftIx += 1
#			(lftRpt,lftOp) = cigarInfo.operations[lftIx]
#
#		# determine how much to absorb from each op based on the op types; and
#		# for each base of the reference consumed, yield the counts
#
#		bothOp = lftOp+rgtOp
#		if (bothOp in ["==","=X","X=","XX"]):
#			rptAbsorbed = min(lftRpt,rgtRpt)
#			for _ in range(rptAbsorbed):
#				if (lftOp == "="): numMatch    -= 1
#				else:              numMismatch -= 1
#				if (rgtOp == "="): numMatch    += 1
#				else:              numMismatch += 1
#				yield (..,numMatch,numMismatch,numIndel)
#			lftRpt -= rptAbsorbed
#			rgtRpt -= rptAbsorbed
#		elif (bothOp == "II"):
#			rptAbsorbed = min(lftRpt,rgtRpt)
#			for _ in range(rptAbsorbed):
#				# (note that numIndel doesn't change)
#				yield (..,numMatch,numMismatch,numIndel)
#			lftRpt -= rptAbsorbed
#			rgtRpt -= rptAbsorbed
#		elif (bothOp == "DD"):
#			# note: no reference consumed, so no yield
#			numIndel -= lftRpt
#			numIndel += rgtRpt
#			lftRpt = rgtRpt = 0
#		elif (bothOp in ["=I","XI"]):
#			rptAbsorbed = min(lftRpt,rgtRpt)
#			for _ in range(rptAbsorbed):
#				if (lftOp == "="): numMatch    -= 1
#				else:              numMismatch -= 1
#				numIndel += 1
#				yield (..,numMatch,numMismatch,numIndel)
#			lftRpt -= rptAbsorbed
#			rgtRpt -= rptAbsorbed
#		elif (bothOp in ["I=","IX"]):
#			rptAbsorbed = min(lftRpt,rgtRpt)
#			for _ in range(rptAbsorbed):
#				numIndel -= 1
#				if (rgtOp == "="): numMatch    += 1
#				else:              numMismatch += 1
#				yield (..,numMatch,numMismatch,numIndel)
#			lftRpt -= rptAbsorbed
#			rgtRpt -= rptAbsorbed
#		elif (bothOp in ["=D","XD"]):
#			# note: no reference consumed, so no yield
#			numIndel += rgtRpt
#			rgtRpt =  0
#		elif (bothOp in ["D=","DX"]):
#			# note: no reference consumed, so no yield
#			numIndel -= lftRpt
#			lftRpt =  0
#		elif (bothOp == "ID"):
#			# note: no reference consumed, so no yield
#			numIndel += rgtRpt
#			rgtRpt =  0
#		elif (bothOp == "DI"):
#			# note: no reference consumed, so no yield
#			numIndel -= lftRpt
#			lftRpt =  0


# split_cigar--
#	Split a cigar string into a list of (count,operation)

class CigarInfo: pass

def split_cigar(cigar):
	if (cigar == "*"): return None

	# split the cigar into a list of (count,operation)

	operations = []
	rpt = []
	for ch in cigar:
		if (ch.isdigit()):
			rpt += [ch]
		else:
			#assert (rpt != []), "bad cigar: \"%s\"" % cigar
			if (rpt == []):
				operations += [(1,ch)]
			else:
				operations += [(int("".join(rpt)),ch)]
			rpt = []
	assert (rpt == []), "bad cigar: \"%s\"" % cigar

	# trim clipping operators from the ends

	startClip = endClip = 0
	if (operations != []):
		(rpt,op) = operations[0]
		if (op == "H"):
			startClip = rpt
			operations = operations[1:]

	if (operations != []):
		(rpt,op) = operations[-1]
		if (op == "H"):
			endClip = rpt
			operations = operations[:-1]

	splitCigar = CigarInfo()
	splitCigar.operations = operations
	splitCigar.startClip  = startClip
	splitCigar.endClip    = endClip
	return splitCigar


# cigar_to_string--
#	Convert a list of (count,operation) into a cigar string
#
# cigarInfo can be a CigarInfo object, or a list of (count,operation) pairs.
# Specifically, if it has an operations attribute, we assume it is a CigarInfo
# object; if not, we assume it is a list

def cigar_to_string(cigarInfo):
	startClip = endClip = 0
	if (hasattr(cigarInfo,"operations")):
		operations = cigarInfo.operations
		startClip  = cigarInfo.startClip
		endClip    = cigarInfo.endClip
	else:
		operations = cigarInfo

	operationsStr = "".join(["%s%s"%((str(rpt) if (rpt>1) else ""),op) for (rpt,op) in operations])
	if (startClip > 0): operationsStr = ("%dH"%startClip) + operationsStr
	if (endClip   > 0): operationsStr = operationsStr + ("%dH"%endClip)

	return operationsStr


# cigar_extent--
#	Determine the extent, the number of bases a cigar covers, in both sequences.
#
# (cigarInfo is as for cigar_to_string)

def cigar_extent(cigarInfo,targetStart=0,queryStart=0,reverseQuery=False):
	if (hasattr(cigarInfo,"operations")):
		operations = cigarInfo.operations
		#startClip  = cigarInfo.startClip  # these don't
		#endClip    = cigarInfo.endClip    # .. affect extent
	else:
		operations = cigarInfo

	qStep = -1 if (reverseQuery) else 1

	(tPos,qPos) = (targetStart,queryStart)
	for (rpt,op) in operations:
		if (op in ["M","X","="]):
			tPos += rpt
			qPos += qStep * rpt
		elif (op == "I"):
			qPos += qStep * rpt
		elif (op == "D"):
			tPos += rpt
		else:
			raise ValueError

	return (tPos,qPos)


# trace_cigar_path--
#	Return the list of (target,query) points in the path represented by a cigar
#	string.
#
# (cigarInfo is as for cigar_to_string)

def trace_cigar_path(cigarInfo,targetStart=0,queryStart=0,reverseQuery=False):
	if (hasattr(cigarInfo,"operations")):
		operations = cigarInfo.operations
		#startClip  = cigarInfo.startClip  # these don't
		#endClip    = cigarInfo.endClip    # .. affect the path
	else:
		operations = cigarInfo

	if (reverseQuery): 
		matchStep  = (1,-1)
		insertStep = (0,-1)
		deleteStep = (-1,0)
	else:
		matchStep  = (1,1)
		insertStep = (0,1)
		deleteStep = (1,0)

	(tPos,qPos) = (targetStart,queryStart)
	path = [(tPos,qPos)]
	for (rpt,op) in operations:
		if (op in ["M","X","="]):
			(tStep,qStep) = matchStep
		elif (op == "I"):
			(tStep,qStep) = insertStep
		elif (op == "D"):
			(tStep,qStep) = deleteStep
		else:
			raise ValueError
		for _ in range(rpt):
			tPos += tStep
			qPos += qStep
			path += [(tPos,qPos)]

	return path


# cigar_before_target--
#	Return (as a list of (count,operation)) the prefix of a cigar string, up to
#	(but not including) a specified target position.
#
# (cigarInfo is as for cigar_to_string)

def cigar_before_target(cigarInfo,targetPos,targetStart=0):
	if (hasattr(cigarInfo,"operations")):
		operations = cigarInfo.operations
		#startClip  = cigarInfo.startClip  # these don't
		#endClip    = cigarInfo.endClip    # .. affect the path
	else:
		operations = cigarInfo

	# follow path until we reach the goal, collecting operations; the final
	# operation is shortened if necessary

	tPos = targetStart
	prefix = []
	for (rpt,op) in operations:
		if (tPos == targetPos):
			break
		assert (tPos < targetPos), "internal error in cigar_before_target, goal passed"
		if (op in ["M","X","=","D"]):
			if (tPos + rpt > targetPos):
				rpt = targetPos - tPos
			prefix += [(rpt,op)]
			tPos += rpt
		elif (op == "I"):
			prefix += [(rpt,op)]
		else:
			raise ValueError

	assert (tPos == targetPos), "internal error in cigar_before_target, goal not reached"

	return prefix


# cigar_after_target--
#	Return (as a list of (count,operation)) the suffix of a cigar string, after
#	(and including) a specified target position.
#
# (cigarInfo is as for cigar_to_string)

def cigar_after_target(cigarInfo,targetPos,targetStart=0):
	if (hasattr(cigarInfo,"operations")):
		operations = cigarInfo.operations
		#startClip  = cigarInfo.startClip  # these don't
		#endClip    = cigarInfo.endClip    # .. affect the path
	else:
		operations = cigarInfo

	# follow path until we reach the goal; the final operation is shortened
	# if necessary, with any leftover part recorded

	tPos = targetStart
	leftover = None
	for (ix,(rpt,op)) in enumerate(operations):
		if (tPos == targetPos):
			break
		assert (tPos < targetPos), "internal error in cigar_after_target, goal passed"
		if (op in ["M","X","=","D"]):
			if (tPos + rpt > targetPos):
				leftover = (tPos+rpt-targetPos,op)
				rpt = targetPos - tPos
			tPos += rpt
		elif (op == "I"):
			pass
		else:
			raise ValueError

	assert (tPos == targetPos), "internal error in cigar_after_target, goal not reached"

	# construct the suffix form the leftover (if any) and the remaining
	# operations

	suffix = [] if (leftover == None) else [leftover]
	suffix += operations[ix:]

	return suffix


# cigar_identity--
#	Compute the nucleotide identity (matches divided by matches+mismatches)
#	represented by a cigar string.
#
# (cigarInfo is as for cigar_to_string)

def cigar_identity(cigarInfo):
	if (hasattr(cigarInfo,"operations")):
		operations = cigarInfo.operations
		#startClip  = cigarInfo.startClip  # these don't
		#endClip    = cigarInfo.endClip    # .. affect the path
	else:
		operations = cigarInfo

	nMatch    = sum([rpt for (rpt,op) in operations if (op in ["M","="])])
	nMismatch = sum([rpt for (rpt,op) in operations if (op == "X")])
	return nMatch / (nMatch+nMismatch)


# cigar_blast_identity--
#	Compute the alignment identity including indels (matches divided by
#	matches+mismatches+indels) represented by a cigar string.
#
# (cigarInfo is as for cigar_to_string)

def cigar_blast_identity(cigarInfo):
	if (hasattr(cigarInfo,"operations")):
		operations = cigarInfo.operations
		#startClip  = cigarInfo.startClip  # these don't
		#endClip    = cigarInfo.endClip    # .. affect the path
	else:
		operations = cigarInfo

	nMatch     = sum([rpt for (rpt,op) in operations if (op in ["M","="])])
	nMutations = sum([rpt for (rpt,op) in operations if (op in ["X","I","D"])])
	return nMatch / (nMatch+nMutations)


# cigar_continuity--
#	Compute the continuity (the fraction of alignment columns that do not
#	contain gaps) represented by a cigar string.
#
# (cigarInfo is as for cigar_to_string)

def cigar_continuity(cigarInfo):
	if (hasattr(cigarInfo,"operations")):
		operations = cigarInfo.operations
		#startClip  = cigarInfo.startClip  # these don't
		#endClip    = cigarInfo.endClip    # .. affect the path
	else:
		operations = cigarInfo

	nMatchOrSub = sum([rpt for (rpt,op) in operations if (op in ["M","=","X"])])
	nColumns    = sum([rpt for (rpt,op) in operations])
	return nMatchOrSub / nColumns


# cigar_nmatch--
#	Count the number of matched bases represented by a cigar string.
#
# (cigarInfo is as for cigar_to_string)

def cigar_nmatch(cigarInfo):
	if (hasattr(cigarInfo,"operations")):
		operations = cigarInfo.operations
		#startClip  = cigarInfo.startClip  # these don't
		#endClip    = cigarInfo.endClip    # .. affect the path
	else:
		operations = cigarInfo

	nMatch = sum([rpt for (rpt,op) in operations if (op in ["M","="])])
	return nMatch

