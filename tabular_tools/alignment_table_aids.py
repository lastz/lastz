#!/usr/bin/env python3
"""
Things built upon the AlignmentTable objects.
"""

# construct_alignment_text--
#	Construct the alignment text, tracing a cigar string's path through the
#	sequences.

def construct_alignment_text(a,cigar,targetLookup,queryLookup):
	tSubstring = targetLookup[a.name1][a.start1:a.end1]
	qSubstring = queryLookup [a.name2][a.start2:a.end2]
	if (a.strand == "-"): qSubstring = reverse_complement(qSubstring)

	text1 = []
	text2 = []

	tIx = qIx = 0
	for (rpt,op) in cigar.operations:
		if (op in ["M","X","="]):
			text1 += tSubstring[tIx:tIx+rpt]
			text2 += qSubstring[qIx:qIx+rpt]
			tIx += rpt
			qIx += rpt
		elif (op == "I"):
			text1 += ["-"] * rpt
			text2 += qSubstring[qIx:qIx+rpt]
			qIx += rpt
		elif (op == "D"):
			text1 += tSubstring[tIx:tIx+rpt]
			text2 += ["-"] * rpt
			tIx += rpt
		else:
			assert (False), "(at line %s) unsupported \"%d%s\" in cigar \"%s\"" \
			              % (("{:,}".format(a.lineNumber)),rpt,op,a.cigarx)

	if (tIx != a.end1 - a.start1):
		assert (False), ("(at line %s) internal error(?);" 
		               + " target %s interval (%d,%d) doesn't match cigar ends (%d,%d)") \
		              % (("{:,}".format(a.lineNumber)),a.name1,a.start1,a.end1,a.start1,a.start1+tIx)

	if (qIx != a.end2 - a.start2):
		assert (False), ("(at line %s) internal error(?);" 
		               + " query %s interval (%d,%d) doesn't match cigar ends (%d,%d)") \
		              % (("{:,}".format(a.lineNumber)),a.name2,a.start2,a.end2,a.start2,a.start2+qIx)

	return ("".join(text1),"".join(text2))


# alignment_score--
#	Compute the score, using lastz default scoring, of an alignment. The
#	alignment is provided as the two aligning strings text1 and text2.
#
# Nota bene: per https://lastz.github.io/lastz/#options_scoring, the first base
#            in a gap incurs the sum of gap open and extend penalties.

# $$$ ideally we'd like to be able to use other scoring schemes, read scores
#     .. from files, etc.

def alignment_score(text1,text2):
	gapOpenPenalty    = 400     # (see note
	gapExtendPenalty  = 30      #  .. above)
	xPenalty          = 1000
	nPenalty          = 100
	substitutionScore = {"AA":   91, "AC":-114, "AG": -31, "AT":-123,
	                     "CA": -114, "CC": 100, "CG":-125, "CT": -31,
	                     "GA":  -31, "GC":-125, "GG": 100, "GT":-114,
	                     "TA": -123, "TC": -31, "TG":-114, "TT":  91}

	# $$$ there's bound to be a faster way to do this

	score = 0
	gapLen1 = gapLen2 = 0
	for (ch1,ch2) in zip(text1.upper(),text2.upper()):
		if (ch1 != "-") and (ch2 != "-"):
			if (gapLen1 > 0):
				score -= gapOpenPenalty + gapLen1*gapExtendPenalty
				gapLen1 = 0
			elif (gapLen2 > 0):
				score -= gapOpenPenalty + gapLen2*gapExtendPenalty
				gapLen2 = 0
			key = ch1+ch2
			if (key in substitutionScore):
				score += substitutionScore[key]
			elif (ch1 == "X") or (ch2 == "X"):
				score -= xPenalty
			else:
				score -= nPenalty
		elif (ch1 == "-"):
			if (gapLen2 > 0):
				score -= gapOpenPenalty + gapLen2*gapExtendPenalty
				gapLen2 = 0
			gapLen1 += 1
		elif (ch2 == "-"):
			if (gapLen1 > 0):
				score -= gapOpenPenalty + gapLen1*gapExtendPenalty
				gapLen1 = 0
			gapLen2 += 1
		else:
			pass  # both nts are a gap -- just ignore it

	if (gapLen1 > 0):
		score -= gapOpenPenalty + gapLen1*gapExtendPenalty
	elif (gapLen2 > 0):
		score -= gapOpenPenalty + gapLen2*gapExtendPenalty

	return score


# reverse_complement--

complementMap = str.maketrans("ACGTSWRYMKBDHVNacgtswrymkbdhvn",
                              "TGCASWYRKMVHDBNtgcaswyrkmvhdbn")

def reverse_complement(nukes):
	return nukes[::-1].translate(complementMap)
