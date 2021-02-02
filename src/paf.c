//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: paf.c
//
//----------
//
// paf--
//	Support for printing alignments in PAF format.
//
// PAF format is a well-established multiple alignment format.  As of Jan/2009,
// a spec for PAF files can be found at
//	http://genome.ucsc.edu/FAQ/FAQformat#format5
//
//----------

//----------
//
// other files
//
//----------

#include <stdlib.h>				// standard C stuff
#define  true  1
#define  false 0
#include <string.h>				// standard C string stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff
#include "diag_hash.h"			// diagonals hashing stuff
#include "identity_dist.h"		// identity distribution stuff
#include "coverage_dist.h"		// query coverage distribution stuff
#include "continuity_dist.h"	// query continuity distribution stuff

#define  paf_owner				// (make this the owner of its globals)
#include "paf.h"				// interface to this module

static int max_digits (s64 num1, s64 num2);

// debugging defines

//#define debugSeq1Beg 858		// if defined, only alignments entirely within
//#define debugSeq1End 1153		// .. this range in sequence 1 are output;  note
								// .. that these positions are origin-zero

//#define snoopPafAlignList		// if this is defined, extra code is added to
								// .. track alignment lists



//----------
//
// print_paf_align_list--
//	Print a list of gapped alignments in paf format.
//
//----------
//
// Arguments:
//	FILE*		f:				The file to print to.
//	alignel*	alignList:		The list of alignments to print.
//	seq*		seq1:			One sequence.
//	seq*		seq2:			Another sequence.
//	int			withComments:	true => print comments as well
//
// Returns:
//	(nothing)
//
//----------

//=== stuff for snoopPafAlignList ===

#ifndef snoopPafAlignList
#define snoopPafAlignList_1 ;
#endif // not snoopPafAlignList

#ifdef snoopPafAlignList

#define snoopPafAlignList_1                                                    \
	fprintf (stderr, "print_paf_align_list  a=%08lX"                           \
	                 "  a->seq1=%08lX  a->seq2=%08lX\n",                       \
	                 (long) a, (long) a->seq1, (long) a->seq2);

#endif // snoopPafAlignList


// print_paf_align_list--

void print_paf_align_list
   (FILE*		f,
	alignel*	alignList,
	seq*		seq1,
	seq*		seq2,
	int			withComments)
	{
	alignel*	a;
	unspos		numer, denom;

	for (a=alignList ; a!=NULL ; a=a->next)
		{
		snoopPafAlignList_1;
		if (withComments)
			{
			unspos height, width, i, j, prevI, prevJ, run;
			u32    opIx;

			// report identity
			alignment_identity (seq1, seq2, a, &numer, &denom);
			fprintf (f, "# identity=" unsposSlashFmt, numer, denom);
			if (denom != 0) fprintf (f, " (%.1f%%)", (100.0*numer) / denom);
			fprintf (f, "\n");

			// report coverage
			alignment_coverage (seq1, seq2, a, &numer, &denom);
			fprintf (f, "# coverage=" unsposSlashFmt, numer, denom);
			if (denom != 0) fprintf (f, " (%.1f%%)", (100.0*numer) / denom);
			fprintf (f, "\n");

			// report continuity
			alignment_continuity (a, &numer, &denom);
			fprintf (f, "# continuity=" unsposSlashFmt, numer, denom);
			if (denom != 0) fprintf (f, " (%.1f%%)", (100.0*numer) / denom);
			fprintf (f, "\n");

			// report alignment path

			fprintf (f, "# cigar=");

			height = a->end1 - a->beg1 + 1;
			width  = a->end2 - a->beg2 + 1;

			opIx = 0;
			for (i=j=0 ; (i< height)||(j<width) ; )
				{
				run = edit_script_run_of_subs (a->script, &opIx);
				if (run > 0)
					{
					fprintf (f, unsposFmt "m", run);
					i += run; j += run;
					}
		
				if ((i < height) || (j < width))
					{
					prevI = i;  prevJ = j;
					edit_script_indel_len (a->script, &opIx, &i, &j);
					if (i > prevI)
						fprintf (f, unsposFmt "d", i - prevI);
					if (j > prevJ)
						fprintf (f, unsposFmt "i", j - prevJ);
					}
				}
			fprintf (f, "\n");
			}

		print_paf_align (f,
		                 seq1, a->beg1-1, a->end1,
		                 seq2, a->beg2-1, a->end2,
		                 a->script, a->s);
		}
	}

//----------
//
// print_paf_align--
//	Print a single gapped alignment in paf format.
//
//----------
//
// Arguments:
//	FILE*		f:				The file to print to.
//	seq*		seq1:			One sequence.
//	unspos		beg1, end1:		Range of positions in sequence 1 (origin 0).
//	seq*		seq2:			Another sequence.
//	unspos		beg2, end2:		Range of positions in sequence 2 (origin 0).
//	editscript*	script:			The script describing the path the alignment
//								.. takes in the DP matrix.
//	score		s:				The alignment's score.
//
// Returns:
//	(nothing)
//
//----------

static char* rcfSuffix[4] = { "", "~", "~", "" };

void print_paf_align
   (FILE*			f,
	seq*			seq1,
	unspos			beg1,
	unspos			end1,
	seq*			seq2,
	unspos			beg2,
	unspos			end2,
	editscript*		script,
	score			s)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part;
	unspos			height, width, i, j, run;
	u32				opIx;
	u8*				p, *q;
	unspos			ix;
	char*			name1, *name2, *pref2, *suff1, *suff2;
	unspos			offset1, offset2, start1, start2;
	unspos			startLoc1, startLoc2;
	unspos			seq1Len, seq2Len, seq1True, seq2True;
	char			strand1, strand2;
	unspos			startI, startJ;
	int				len1, len2, nameW, startW, endW, lenW;

#ifdef debugSeq1Beg
	if ((beg1 < debugSeq1Beg) || (end1 > debugSeq1End)) return;
#endif // debugSeq1Beg

	beg1++; // (internally, we want origin 1, inclusive)
	beg2++;

	height = end1 - beg1 + 1;
	width  = end2 - beg2 + 1;

	// report diagonal

	if (paf_dbgReportDiag)
		fprintf (f, "# diagonal=" sgnposFmt "\n", diagNumber(beg1,beg2));

	//////////
	// figure out position offsets and names
	//////////

	if (sp1->p == NULL)		// sequence 1 is not partitioned
		{
		name1 = (seq1->useFullNames)? seq1->header : seq1->shortHeader;
		if ((name1 == NULL) || (name1[0] == 0)) name1 = "seq1";
		offset1   = 0;
		startLoc1 = seq1->startLoc;
		seq1Len   = seq1->len;
		seq1True  = seq1->trueLen;
		}
	else					// sequence 1 is partitioned
	 	{
		part = lookup_partition (seq1, beg1-1);
		name1     = &sp1->pool[part->header];
		offset1   = part->sepBefore + 1;
		startLoc1 = part->startLoc;
		seq1Len   = part->sepAfter - offset1;
		seq1True  = part->trueLen;
		}

	if (sp2->p == NULL)		// sequence 2 is not partitioned
		{
		name2 = (seq2->useFullNames)? seq2->header : seq2->shortHeader;
		if ((name2 == NULL) || (name2[0] == 0)) name2 = "seq2";
		offset2   = 0;
		startLoc2 = seq2->startLoc;
		seq2Len   = seq2->len;
		seq2True  = seq2->trueLen;
		}
	else					// sequence 2 is partitioned
	 	{
		part = lookup_partition (seq2, beg2-1);
		name2     = &sp2->pool[part->header];
		startLoc2 = part->startLoc;
		offset2   = part->sepBefore + 1;
		seq2Len   = part->sepAfter - offset2;
		seq2True  = part->trueLen;
		}

	//////////
	// print summary line
	//////////

	fprintf (f, "a score=" scoreFmt "\n", s);

	//////////
	// print aligning path in sequence 1
	//////////

	// figure out fields and widths

	pref2 = ((paf_distinguishNames) && (strcmp (name1, name2) == 0))? "~" : "";
	suff1 = rcfSuffix[seq1->revCompFlags];
	suff2 = rcfSuffix[seq2->revCompFlags];

	if ((seq1->revCompFlags & rcf_rev) == 0)
		{
		start1  = beg1-1 - offset1 + startLoc1;
		strand1 = '+';
		}
	else
		{
		start1  = beg1-1 - offset1 + seq1True+2 - (startLoc1 + seq1Len);
		strand1 = '-';
		}
	if ((seq2->revCompFlags & rcf_rev) == 0)
		{
		start2  = beg2-1 - offset2 + startLoc2;
		strand2 = '+';
		}
	else
		{
		start2  = beg2-1 - offset2 + seq2True+2 - (startLoc2 + seq2Len);
		strand2 = '-';
		}

	len1  =                  strlen (name1) + strlen (suff1);
	len2  = strlen (pref2) + strlen (name2) + strlen (suff2);
	nameW = (len1 >= len2)? len1 : len2;

	startW = max_digits (start1, start2);
	endW   = max_digits (end1+1-beg1, end2+1-beg2);
	lenW   = max_digits (seq1True, seq2True);

	// print aligning path in sequence 1 (non-printables are printed as '*'
	// but such should never be seen unless there is a problem elsewhere)

	fprintf (f, "s %s%s%*s" unsposStarFmt " " unsposStarFmt " %c " unsposStarFmt " ",
	            name1, suff1, nameW+1-len1, " ",
	            startW, start1-1, endW, end1+1-beg1, strand1, lenW, seq1True);

	opIx = 0;
	for (i=j=0 ; (i<height)||(j<width) ; )
		{
		// handle the next run

		run = edit_script_run_of_subs (script, &opIx);

		p = seq1->v+beg1+i-1;
		q = seq2->v+beg2+j-1;
		for (ix=0 ; ix<run ; ix++)
			{ fprintf (f, "%c", dna_toprint(*p));  p++;  q++; }

		i += run; j += run;

		// handle the next indel

		if ((i < height) || (j < width))
			{
			startI = i;  p = seq1->v+beg1+i-1;
			startJ = j;  q = seq2->v+beg2+j-1;

			edit_script_indel_len (script, &opIx, &i, &j);

			if (i != startI)
				{
				for ( ; startI<i ; startI++)
					{ fprintf (f, "%c", dna_toprint(*p));  p++; }
				}

			if (j != startJ)
				{
				for ( ; startJ<j ; startJ++)
					{ fprintf (f, "-");  q++; }
				}
			}
		}

	fprintf (f, "\n");

	//////////
	// print aligning path in sequence 2
	//////////

	fprintf (f, "s %s%s%s%*s" unsposStarFmt " " unsposStarFmt " %c " unsposStarFmt " ",
	            pref2, name2, suff2, nameW+1-len2, " ",
	            startW, start2-1, endW, end2+1-beg2, strand2, lenW, seq2True);

	opIx = 0;
	for (i=j=0 ; (i<height)||(j<width) ; )
		{
		// handle the next run

		run = edit_script_run_of_subs (script, &opIx);

		p = seq1->v+beg1+i-1;
		q = seq2->v+beg2+j-1;
		for (ix=0 ; ix<run ; ix++)
			{ fprintf (f, "%c", dna_toprint(*q));  p++;  q++; }

		i += run; j += run;

		// handle the next indel

		if ((i < height) || (j < width))
			{
			startI = i;  p = seq1->v+beg1+i-1;
			startJ = j;  q = seq2->v+beg2+j-1;

			edit_script_indel_len (script, &opIx, &i, &j);

			if (i != startI)
				{
				for ( ; startI<i ; startI++)
					{ fprintf (f, "-");  p++; }
				}

			if (j != startJ)
				{
				for ( ; startJ<j ; startJ++)
					{ fprintf (f, "%c", dna_toprint(*q));  q++; }
				}
			}
		}

	fprintf (f, "\n\n");
	}

//----------
//
// max_digits--
//	Figure out the number of digits required to print either of two numbers.
//
//----------
//
// Arguments:
//	s64 num1,num2: The two numbers.
//
// Returns:
//	The number of characters required for printing.
//
//----------

static int max_digits
   (s64		num1,
	s64		num2)
	{
	int w1 = snprintf (NULL, 0, s64Fmt, num1);
	int w2 = snprintf (NULL, 0, s64Fmt, num2);
	return (w1 >= w2)? w1 : w2;
	}
