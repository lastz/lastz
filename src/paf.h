//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: paf.h
//
//----------

#ifndef paf_H					// (prevent multiple inclusion)
#define paf_H

// other files

#include <stdio.h>				// standard C i/o stuff
#include <stdarg.h>				// standard C variable argument list stuff
#include "utilities.h"			// utility stuff
#include "sequences.h"			// sequence stuff
#include "edit_script.h"		// alignment edit script stuff

// establish ownership of global variables

#ifdef paf_owner
#define global
#else
#define global extern
#endif

// "deep link" control variable access

#ifdef paf_owner
int paf_distinguishNames = false;	// true => add a "~" prefix to the second
									// sequence name when names are identical
int paf_dbgReportDiag = false;		// true => report diagonal as a paf comment
#else
global int paf_distinguishNames;
global int paf_dbgReportDiag;
#endif

//----------
//
// prototypes for routines in paf.c
//
//----------

void print_paf_job_header (FILE* f,
                           char* programName, char* args, scoreset* scoring,
                           sthresh* hspThreshold, sthresh* gappedThreshold,
                           score xDrop, score yDrop,
                           int withComments);
void print_paf_job_footer (FILE* f);
void print_paf_header     (FILE* f, seq* seq1, seq* seq2);
void print_paf_align_list (FILE* f, alignel* alignList, seq* seq1, seq* seq2,
                           int withComments);
void print_paf_align      (FILE* f,
                           seq* seq1, unspos beg1, unspos end1,
                           seq* seq2, unspos beg2, unspos end2,
                           editscript* script, score s);
void print_paf_match      (FILE* f,
                           seq* seq1, unspos pos1,
                           seq* seq2, unspos pos2, unspos length,
                           score s, int withComments);
void print_paf_comment    (FILE* f, const char* format, ...);
void vprint_paf_comment   (FILE* f, const char* format, va_list args);

#undef global
#endif // paf_H
