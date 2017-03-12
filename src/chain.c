//-------+---------+---------+---------+---------+---------+---------+--------=
//
// File: chain.c
//
//----------
//
// chain--
//	Find the highest scoring chain in a set of gap-free alignments.  Each
//	segment in the chain will begin strictly before the start of the next
//	segment.  This is (expected to be) the most parsimonious subset of the
//	gap-free alignments, assuming there actual orthology contains no
//	inversions.
//
// The algorithm finds, for each segment, the highest scoring chain that ends
// with that segment.  Segments are scanned in an order (by increasing start
// in sequence 1) that guarantees all possible predecessor chains have been
// found and scored before that segment is considered.  Upon completion, the
// chain is recovered by backtracking from its end segment.
//
// A chain's score is the sum of its segment scores minus the sum of penalties
// for the gaps between segments.  The caller must provide a function to compute
// those penalties.  See note (1) of reduce_to_chain() for more details.
//
// To facilitate the search for valid predecessors, a K-d tree is used.  See
// the header of build_kd_tree() for more details on the tree implementation.
//
// References:
//
//	[1] Multidimensional Binary Search Trees Used for Associative Searching.
//	    Jon Louis Bentley, Commun. ACM 18(9): 509-517 (1975).
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
#include <stdio.h>				// standard C i/o stuff
#include <string.h>				// standard C string stuff
#include "build_options.h"		// build options
#include "utilities.h"			// utility stuff
#include "dna_utilities.h"		// dna/scoring stuff
#include "sequences.h"			// sequence stuff

#define  chain_owner			// (make this the owner of its globals)
#include "chain.h"				// interface to this module

// debugging defines

//#define snoopBatches			// if this is defined, extra code is added to
								// .. track the chaining of HSPs in per-
								// .. partition batches

//----------
//
// private global data
//
//----------

typedef double bigscore;

typedef struct kdinfo
	{
	score		diagPen, antiPen;		// chain gap penalties
	int			scale;					// chain score scale factor
	bigscore*	chainScore;				// array of the best score for any
										// .. chain ending with a given segment
	u32*		perm, *invPerm;			// permuation of segments (and inverse)
	segment*	seg;					// array of segments
	segment*	query;					// query segment
	unspos		x, y;					// query's cartesian point
	sgnpos		diag;					// .. x = pos1;  y = pos2;  diag = x-y
	chainer		connect;				// chain connection penalty function
	} kdinfo;


typedef struct bestpred
	{
	u32			num;					// index of a predecessor segment;  the
										// .. value noPred indicates there is no
										// .. predecessor
	bigscore	contrib;				// score of the chain ending at that
										// .. segment, inluding penalty to
										// .. connect it to the query segment
	} bestpred;

#define noPred ((u32) -1)


typedef struct kdnode
	{
	int			isBucket;				// true => this node is a bucket/leaf
	u32			loIx, hiIx;				// isBucket is true
										// ..   index range of the segments in
										// ..   this leaf
										// isBucket is false
										// ..   hiIx is index corresponding to
										// ..   cutVal

	// the following fields are onlyl valid if isBucket is false

	sgnpos		cutVal;					// value (along appropriate axis) which
										// .. separates lower and upper children
	bigscore	maxChainScore;			// the highest score for any chain
										// .. ending at a segment in this
										// .. subtree
	struct kdnode* loSon, *hiSon;		// pointers to child nodes
	} kdnode;

#define bucketSize 3	// max number of entries we'll place in a bucket node

#define valid_kdnode(p) \
	(((p)->isBucket) || (((p)->loSon != NULL) && ((p)->hiSon != NULL)))

#define perm_swap(kdi,p,q) \
  { u32 t; t = kdi->perm[p]; kdi->perm[p] = kdi->perm[q]; kdi->perm[q] = t; }

// projection-- figure out spatial position of segment i along the current axis

#define projection(i,axis,kdi) \
	((axis == 0) ? (((sgnpos)kdi->seg[kdi->perm[i]].pos1) - ((sgnpos)kdi->seg[kdi->perm[i]].pos2)) \
	             :  ((sgnpos)kdi->seg[kdi->perm[i]].pos2))


// segment batches-- partitioning of a segment table

typedef struct segbatch
	{
	u32		start;				// index (into a segment table) of the first
								// .. entry in a batch
	u32		end;				// index (into a segment table) of the first
								// .. entry NOT in a batch (i.e. the one after
								// .. the last entry).
	partition*	part1;			// sequence partitions that "contain" this
	partition*	part2;			// .. batch;  either of these can be NULL if we
								// .. aren't dealing with partitions in the
								// .. corresponding sequence
	} segbatch;

typedef struct sbtable
	{
	u32		size;				// the number of entries allocated for batch[]
	u32		len;				// the number of batches (the number of entries
								// .. actually used)
	segbatch batch[1];			// the batch table (variable-length array)
	} sbtable;

#define sbtable_bytes(size) (sizeof(sbtable) + (((size)-1)*sizeof(segbatch)))

//----------
//
// prototypes for private functions
//
//----------

static kdnode*  build_kd_tree       (u32 lo, u32 hi, int axis,
                                     const kdinfo* const kdi);
static void     free_kd_tree        (kdnode* subtree);
static void     dump_kd_tree        (FILE* f, int points, kdnode* root,
                                     const kdinfo* const kdi);
static u32      partition_segments  (u32 lo, u32 hi, int axis,
                                     const kdinfo* const kdi);
static bestpred best_predecessor    (kdnode* subtree, int axis,
                                     bigscore lowerBound,
                                     bestpred bp, const kdinfo* const kdi);
static void     propagate_max_score (kdnode* subtree, bigscore s, u32 ix);

//----------
//
// try_reduce_to_chain--
//	This routine is a wrapper for reduce_to_chain(), to handle cases that arise
//	when either of sequence 1 or sequence 2 is partitioned.
//
//----------
//
// Arguments:
//	seq*		seq1:		The first sequence(s).
//	seq*		seq2:		The second sequence(s).
//	(the rest are the same as for reduce_to_chain)
//
// Returns:
//	The score of the best chain, unscaled;  zero if there's some problem.
//	
//
//----------

//=== stuff for snoopBatches ===

#ifndef snoopBatches
#define debugSnoopBatches_1 ;
#endif // not snoopBatches

#ifdef snoopBatches

#define debugSnoopBatches_1                                                   \
	{                                                                         \
	partition* batPart1 = chainBatches->batch[batIx].part1;                   \
	partition* batPart2 = chainBatches->batch[batIx].part2;                   \
	fprintf (stderr, "batch[%u] %u..%u",                                      \
					 batIx, startSegIx, endSegIx-1);                          \
	if (batPart1 != NULL)                                                     \
		fprintf (stderr, " seq1: " unsposFmt ".." unsposFmt,                  \
						 batPart1->sepBefore+1, batPart1->sepAfter);          \
	if (batPart2 != NULL)                                                     \
		fprintf (stderr, " seq2: " unsposFmt ".." unsposFmt,                  \
						 batPart2->sepBefore+1, batPart2->sepAfter);          \
	fprintf (stderr, "\n");                                                   \
	}
#endif // snoopBatches


// private data for try_reduce_to_chain

static sbtable* chainBatches = NULL;


// try_reduce_to_chain--

score try_reduce_to_chain
   (seq*			seq1,
	seq*			seq2,
	segtable*		st,
	score			diagPen,
	score			antiPen,
	int				scale,
	chainer			connect)
	{
	seqpartition*	sp1 = &seq1->partition;
	seqpartition*	sp2 = &seq2->partition;
	partition*		part1, *part2;
	segment*		seg, *seg2;
	segtable		stSubset;
	u32				entriesNeeded;
	size_t			bytesNeeded;
	u32				partIx1, partIx2, segIx, segIx2;
	u32				batIx, startSegIx, endSegIx, segsKept;
	unspos			pEnd1, pEnd2;
	bigscore		chainScore, best;

	// if neither sequence is partitioned, just pass the segments along to
	// reduce_to_chain()

	if ((sp1->p == NULL)	// sequence 1 is not partitioned
	 && (sp2->p == NULL))	// sequence 2 is not partitioned
		return reduce_to_chain (st,diagPen,antiPen,scale,connect);

	// allocate the batch table, or resize it if it's not big enough

	entriesNeeded = 1;
	if (sp1->p != NULL) entriesNeeded *= sp1->len;	// sequence 1 is partitioned
	if (sp2->p != NULL) entriesNeeded *= sp2->len;	// sequence 2 is partitioned
	bytesNeeded = sbtable_bytes (entriesNeeded);

	if (chainBatches == NULL)
		{
		chainBatches = (sbtable*) malloc_or_die ("try_reduce_to_chain", bytesNeeded);
		chainBatches->size = entriesNeeded;
		}
	else if (entriesNeeded > chainBatches->size)
		{
		chainBatches = (sbtable*) realloc_or_die ("try_reduce_to_chain", chainBatches, bytesNeeded);
		chainBatches->size = entriesNeeded;
		}

	// loop over all partitions in both target and query, collecting batches of
	// segments that are confined to a single partition;  unfortunately each of
	// the three cases (regarding target and/or query partitioned) has to be
	// handled separately

	if ((sp1->p != NULL)	// sequence 1 is partitioned
	 && (sp2->p == NULL))	// sequence 2 is not partitioned
		{
		sort_segments (st, qSegmentsByPos1);

		batIx = 0;
		segIx = 0;  seg = &st->seg[segIx];
		for (partIx1=0 ; partIx1<sp1->len ; partIx1++)
			{
			if (segIx >= st->len) break;

			part1 = &sp1->p[partIx1];
			pEnd1 = part1->sepAfter;
			if (pEnd1 < seg->pos1 + seg->length) continue;

			chainBatches->batch[batIx].part1 = part1;
			chainBatches->batch[batIx].part2 = NULL;
			chainBatches->batch[batIx].start = segIx++;  seg++;
			while ((segIx < st->len) && (pEnd1 >= seg->pos1 + seg->length))
				{ segIx++;  seg++; }

			chainBatches->batch[batIx].end = segIx;
			batIx++;
			}
		chainBatches->len = batIx;
		}

	else if ((sp1->p == NULL)	// sequence 1 is not partitioned
	      && (sp2->p != NULL))	// sequence 2 is partitioned
		{
		sort_segments (st, qSegmentsByPos2);

		batIx = 0;
		segIx = 0;  seg = &st->seg[segIx];
		for (partIx2=0 ; partIx2<sp2->len ; partIx2++)
			{
			if (segIx >= st->len) break;

			part2 = &sp2->p[partIx2];
			pEnd2 = part2->sepAfter;
			if (pEnd2 < seg->pos2 + seg->length) continue;

			chainBatches->batch[batIx].part1 = NULL;
			chainBatches->batch[batIx].part2 = part2;
			chainBatches->batch[batIx].start = segIx++;  seg++;
			while ((segIx < st->len) && (pEnd2 >= seg->pos2 + seg->length))
				{ segIx++;  seg++; }

			chainBatches->batch[batIx].end = segIx;
			batIx++;
			}
		chainBatches->len = batIx;
		}

	else if ((sp1->p != NULL)	// sequence 1 is partitioned
	      && (sp2->p != NULL))	// sequence 2 is partitioned
		{
		sort_segments (st, qSegmentsByPos1);

		batIx = 0;
		segIx = 0;  seg = &st->seg[segIx];
		for (partIx1=0 ; partIx1<sp1->len ; partIx1++)
			{
			if (segIx >= st->len) break;

			part1 = &sp1->p[partIx1];
			pEnd1 = part1->sepAfter;
			if (pEnd1 < seg->pos1 + seg->length) continue;

			startSegIx = segIx++;  seg++;
			while ((segIx < st->len) && (pEnd1 >= seg->pos1 + seg->length))
				{ segIx++;  seg++; }
			endSegIx = segIx;

			sort_some_segments (st, startSegIx, endSegIx, qSegmentsByPos2);

			segIx2 = startSegIx;  seg2 = &st->seg[segIx2];
			for (partIx2=0 ; partIx2<sp2->len ; partIx2++)
				{
				if (segIx2 >= endSegIx) break;
	
				part2 = &sp2->p[partIx2];
				pEnd2 = part2->sepAfter;
				if (pEnd2 < seg2->pos2 + seg2->length) continue;
	
				chainBatches->batch[batIx].part1 = part1;
				chainBatches->batch[batIx].part2 = part2;
				chainBatches->batch[batIx].start = segIx2++;  seg2++;
				while ((segIx2 < endSegIx) && (pEnd2 >= seg2->pos2 + seg2->length))
					{ segIx2++;  seg2++; }
	
				chainBatches->batch[batIx].end = segIx2;
				batIx++;
				}
			}
		chainBatches->len = batIx;
		}

	// perform chaining over each batch of segments, treating each as a
	// separate chaining problem

	best = 0;

	for (batIx=0 ; batIx<chainBatches->len ; batIx++)
		{
		// create a subset of the segment table, corresponding to this batch

		startSegIx = chainBatches->batch[batIx].start;
		endSegIx   = chainBatches->batch[batIx].end;
		debugSnoopBatches_1

		subset_segment_table (st, startSegIx, endSegIx, &stSubset);

		// chain the subset

		chainScore = reduce_to_chain (&stSubset, diagPen, antiPen, scale, connect);
		if (chainScore > best) best = chainScore;

		// mark the segments that are part of the chain as to-be-kept, and the
		// others as to-be-filtered;  note that reduce_to_chain has brought the
		// chained segments to the front of the subset -- the rest are
		// essentially garbage (and don't necessarily represent the excluded
		// segments of the subset)

		segsKept = stSubset.len;
		for (segIx=startSegIx ; segIx<startSegIx+segsKept ; segIx++)
			st->seg[segIx].filter = false;
		for ( ; segIx<endSegIx ; segIx++)
			st->seg[segIx].filter = true;
		}

	// perform the final filtering step, unaware of batching

	filter_marked_segments (st);

	return best;
	}

//----------
//
// reduce_to_chain--
//	Find the highest scoring chain, in which each segment in the chain begins
//	strictly before the start of the next segment.
//
// WARNING: External modules should usually NOT call this function directly,
//          but should instead call try_reduce_to_chain.  The exception is that
//          if the caller can guarantee that all segments are within the same
//          partition (in both target and query), they can call this directly.
//
// A chain is a series of segments, where each segment in the chain (other than
// the last), begins strictly before the start of the next.  A chain's score is
// scale times the sum of segment scores minus the sum of penalties for the gaps
// between segments:
//		connect (segment_i, segment_(i+1), scale)
// the last sum is taken over all segments in the chain except the last).
//
//----------
//
// Arguments:
//	segtable*	st:			The segments on which to operate.
//	score		diagPen:	Chaining penalty;  see notes (1) and (3).
//	score		antiPen:	Chaining penalty;  see notes (1) and (3).
//	int			scale:		Scaling constant;  see note (2).
//	chainer		connect:	Chain connection penalty function;  see note above,
//							.. and description of arguments in chain.h
//
// Returns:
//	The score of the best chain, unscaled;  zero if there's some problem.  Note
//	that the segment table is modified in place, with segments belonging to the
//	chain brought to the front, and the table shortened by modifying st->len.
//
//----------
//
// Notes:
//	(1)	The parameters diagPen and antiPen permit us to deduce useful
//		inequalities about chain scores.  Namely, let segment_i and segment_j
//		be segments on diagonals diag_i and diag_j, and set
//			diff = diag_j - diag_i
//		Then diagPen and antiPen are required to satisfy:
//			if diff >= 0, then connect(segment_i,segment_j) >= diff*diagPen
//		and
//			if diff < 0, then connect(segment_i,segment_j) >= -diff*antiPen
//
//	(2)	In effect, scale permits integer arithmetic to be used with very small
//		gap penalties, since the computed chain also maximizes the sum of the
//		segment scores minus the sum of
//			connect(segment_i, segment_(i+1), scale)/scale.
//
//	(3) diagPen and antiPen are considered to have already been scaled.  We
//		only apply scale to the segment substitution scores.
//
//----------

#define debugChaining_1                                                      \
	if (chain_dbgChaining)                                                   \
		fprintf (stderr,                                                     \
		         "chaining [%d] " unsposSlashFmt                             \
		         "\tdiag=" sgnposFmt                                         \
		         "\tsegscore=" scoreFmtSimple "\n",                          \
				 i, kdi.x, kdi.y, kdi.diag, kdi.query->s);

#define debugChaining_2                                                      \
	if (chain_dbgChaining)                                                   \
		{                                                                    \
		if (bp.num == noPred)                                                \
			fprintf (stderr, "  pred=(none)\n");                             \
		else                                                                 \
			fprintf (stderr, "  pred=%u query=%.2f contrib=%.2f score=%.2f\n", \
			                 bp.num, queryContrib, bp.contrib,               \
			                 kdi.chainScore[i]);                             \
		}

#define debugChaining_3                                                      \
	if (chain_dbgChaining)                                                   \
		{                                                                    \
		if (bestEnd == noPred)                                               \
			fprintf (stderr, "best=(none)\n");                               \
		else                                                                 \
			fprintf (stderr, "best=%u score=%.2f\n", bestEnd, best);         \
		}


score reduce_to_chain
   (segtable*	st,
	score		diagPen,
	score		antiPen,
	int			scale,
	chainer		connect)
	{
	kdinfo		kdi;
	kdnode*		root;
	u32*		chain;
	bigscore	best, queryContrib;
	u32			bestEnd;
	u32			i, n;
	bestpred	bp;
	segment*	p;

	if (st == NULL) return 0;

	n = st->len;
	if (n == 0) return 0;

	chain_add_stat (numAnchors, n);

	// sort segments by pos1, so that the predecessor search loop is guaranteed
	// to score all possible predecessors of any segment before it considers
	// that segment

	sort_segments (st, qSegmentsByPos1);

	// initialize 'global' data

	kdi.connect    = connect;
	kdi.seg        = st->seg;
	kdi.perm       = malloc_or_die ("reduce_to_chain perm",       n*sizeof(u32));
	kdi.invPerm    = malloc_or_die ("reduce_to_chain invPerm",    n*sizeof(u32));
	kdi.chainScore = zalloc_or_die ("reduce_to_chain chainScore", n*sizeof(bigscore));
	kdi.diagPen    = diagPen;
	kdi.antiPen    = antiPen;
	kdi.scale      = scale;

	// build the K-d tree;  as part of this process the segments are permuted
	// (by use of the perm[] array), and we compute the inverse of that
	// permutation to aid later access to the segments

	for (i=0 ; i<n ; i++)					// build the identity permutation,
		kdi.perm[i] = i;					// .. mapping node numbers to
											// .. segments

	root = build_kd_tree (0, n-1, 1, &kdi);	// build the K-d tree (alters
											// .. kdi.perm)

	for (i=0 ; i<n ; i++)					// compute the inverse permutation
		 kdi.invPerm[kdi.perm[i]] = i;

	if (chain_dbgDumpTree)
		dump_kd_tree (stderr, n, root, &kdi);

	// for each segment, find the best chain ending at that segment;  the array
	// chain[] provides the path from any segment back through its chain;  any
	// segment i for which best_predecessor() finds no positive scoring
	// predecessor (such as those that have no predecesssor at all) will have
	// chain[i] = noPred and terminate backtracking

	chain = malloc_or_die ("reduce_to_chain chain", n*sizeof(u32));

	best    = 0;
	bestEnd = noPred;
	for (i=0 ; i<n ; i++)
		{
		kdi.query = &kdi.seg[i];
		kdi.x     = kdi.query->pos1;
		kdi.y     = kdi.query->pos2;
		kdi.diag  = ((sgnpos) kdi.x) - ((sgnpos) kdi.y);
		debugChaining_1;

		bp.num     = noPred;
		bp.contrib = 0;
		bp = best_predecessor (root, 1, 0, bp, &kdi);
		queryContrib = ((bigscore) kdi.query->s) * ((bigscore) kdi.scale);
		kdi.chainScore[i] = queryContrib + bp.contrib;
		debugChaining_2;

		if (kdi.chainScore[i] > best)
			{ best = kdi.chainScore[i];  bestEnd = i; }
		chain[i] = bp.num;
		propagate_max_score (root, kdi.chainScore[i], kdi.invPerm[i]);
		}

	debugChaining_3;

	// get rid of non-chain segments

	for (p=st->seg ; ((u32)(p-st->seg))<st->len ; p++)
		p->filter = true;

	for (i=bestEnd ; i!=noPred ; i = chain[i])
		(kdi.seg+i)->filter = false;

	filter_marked_segments (st);
	chain_add_stat (numSegments, st->len);

	// scale back best score

	if (dna_utilities_scoreType == 'I')
		{
		best = (best / scale) + 0.5;	// best /= scale, rounded off
		if (best > bestPossibleScore)	// .. and clipped
			best = bestPossibleScore;
		}
	else
		best /= scale;

	free_if_valid ("reduce_to_chain perm",       kdi.perm);
	free_if_valid ("reduce_to_chain invPerm",    kdi.invPerm);
	free_if_valid ("reduce_to_chain chainScore", kdi.chainScore);
	free_if_valid ("reduce_to_chain chain",      chain);
	free_kd_tree  (root);

	return best;
	}

//----------
//
// build_kd_tree--
//	Build segments into a K-d tree.
//
// Standard K-d tree implimentation (for K=2), such as might be found in
// reference [1].  The points are partitioned into two sets, split by the
// a value along one axis.  Each of those sets is in turn split again, along
// the other axis, and so on, until all sets are small enough.  "Small enough"
// is defined by bucketSize.  The two dimensional axes are y (sequence pos2)
// and diagonal (sequence pos1-pos2).
//
//----------
//
// Arguments:
//	u32		lo,hi:	range of entries (of kdi->seg[], indexed by kdi->perm[]) to
//					.. build a tree of;  these are inclusive (i.e. there are
//					.. hi+1-lo entries)
//	int		axis:	which dimension/axis to partition (at the top level)
//					  0 => diagonal (pos1 - pos2)
//					  1 => pos2
//	kdinfo* kdi:	'Global' control variables.
//
// Returns:
//	The root of the tree.  kdi->perm[] is modified so that entries in
//	kdi->seg[kdi->perm[]] agree with the tree.
//
//----------

static kdnode* build_kd_tree
   (u32		lo,
	u32		hi,
	int		axis,
	const kdinfo* const kdi)
	{
	kdnode*	p;
	u32		m;

	p = zalloc_or_die ("build_kd_tree", sizeof(kdnode));
	p->maxChainScore = 0;

	if (hi+1-lo <= bucketSize)		// the range is small enough to fit in one
		{							// .. node
		p->isBucket = true;
		p->loIx = lo;
		p->hiIx = hi;
		}
	else							// the range is two big for one node, split
		{							// .. it into two subtrees
		p->isBucket = false;
		m = partition_segments (lo, hi, axis, kdi);
		p->cutVal = projection (m, axis, kdi);
		p->hiIx = m;
		p->loSon = build_kd_tree (lo,  m,  1-axis, kdi);
		p->hiSon = build_kd_tree (m+1, hi, 1-axis, kdi);
		}

	if (p == NULL)
		suicide ("(in build_kd_tree, p == NULL)");
	if (!valid_kdnode(p))
		suicide ("(in build_kd_tree, p is not a valid kdnode)");

	return p;
	}

//----------
//
// free_kd_tree--
//	Dispose of the memory allocated for a K-d tree.
//
//----------
//
// Arguments:
//	kdnode*		subtree:	The K-d (sub)tree to dispose of.
//
// Returns:
//	(nothing)
//
//----------

// $$$ we could eliminate tail recursion here

static void free_kd_tree
   (kdnode* subtree)
	{
	if (subtree->isBucket)
		free_if_valid ("free_kd_tree leaf", subtree);
	else
		{
		free_kd_tree (subtree->loSon);
		free_kd_tree (subtree->hiSon);
		free_if_valid ("free_kd_tree node", subtree);
		}
	}

//----------
//
// dump_kd_tree--
//	Dump a K-d tree to a file (for debugging).
//
//----------
//
// Arguments:
//	FILE*	f:		The file to print to.
//	int		points:	The number of points in the tree.
//	kdnode*	root:	The K-d tree to print.
//	kdinfo* kdi:	'Global' control variables.
//
// Returns:
//	(nothing)
//
//----------

static void dump_kd_subtree (int indent, kdnode* subtree, int axis);

static FILE*         dksFile;
static const kdinfo* dksKdi;
static int           dksIndexWidth;

static void dump_kd_tree
   (FILE*	f,
	int		points,
	kdnode* root,
	const kdinfo* const kdi)
	{
	int		i;

	for (dksIndexWidth=1,i=points-1 ; i>9 ; i/=10) dksIndexWidth++;

	dksFile = f;
	dksKdi  = kdi;
	dump_kd_subtree (0, root, 1);
	}

static void dump_kd_subtree
   (int			indent,
	kdnode* 	subtree,
	int			axis)
	{
	u32			i, j;
	segment*	s;

	if (!subtree->isBucket)
		{
		dump_kd_subtree (indent+2, subtree->loSon, 1-axis);
		fprintf (dksFile, " %*s %*s%s<=" sgnposFmt "\n", indent, "",
		                  dksIndexWidth, "",
		                  (axis==0)? "x-y" : "x", subtree->cutVal);
		dump_kd_subtree (indent+2, subtree->hiSon, 1-axis);
		return;
		}

	for (i=subtree->loIx ; i<=subtree->hiIx ; i++)
		{
		j = dksKdi->perm[i];
		s = &dksKdi->seg[j];
		fprintf (dksFile, "[%*d]%*s" unsposSlashFmt "\n",
		                  dksIndexWidth, i, indent, "", s->pos1, s->pos2);
		}

	}

//----------
//
// partition_segments--
//	Partition a list of segments into two sets, split by a value along one axis.
//
// A 'pivot' value is desginated as the median (along the specified axis) of
// the first and last segment, and the one at the middle of the list.  Segments
// below the pivot are moved to the first part of the list, and segments above
// it are moved to the last part, with the pivot between them.
//
//----------
//
// Arguments:
//	u32		lo,hi:	range of entries (of kdi->seg[], indexed by kdi->perm[]) to
//					.. build a tree of;  these are inclusive (i.e. there are
//					.. hi+1-lo entries)
//	int		axis:	dimension/axis to partition on
//					  0 => diagonal (pos1 - pos2)
//					  1 => pos2
//	kdinfo* kdi:	'Global' control variables.
//
// Returns:
//	The index of the pivot (m).  kdi->perm[] is modified so that entries in
//	kdi->seg[kdi->perm[]] satsify lo..m-1 <= m <= m+1..hi.
//
//----------

static u32 partition_segments
   (u32		lo,
	u32		hi,
	int		axis,
	const kdinfo* const kdi)
	{
	u32		m, i, j;
	sgnpos	a, b, c, pivot;

	if (hi - lo < 2)
		suicidef ("partition: cannot happen (" unsposCommaFmt ")", lo, hi);

	while (true)
		{
		// find the pivot and move it to the front;  we use the median of the
		// lower, middle, and upper values as the pivot

		m = (lo+hi)/2;
		a = projection (lo, axis, kdi);
		b = projection (m,  axis, kdi);
		c = projection (hi, axis, kdi);
	
		if (((a <= b) && (b <= c)) || ((c <= b) && (b <= a)))
			{ perm_swap (kdi,lo,m);  pivot = b; }
		else if (((a <= c) && (c <= b)) || ((b <= c) && (c <= a)))
			{ perm_swap (kdi,lo,hi);  pivot = c; }
		else
			pivot = a;
	
		// move smaller entries to front, larger to back
	
		i = lo;
		j = hi+1;
		while (i < j)
			{
			// search forward for a large entry
			for (i++ ; (i<=hi)&&(projection(i,axis,kdi)<=pivot) ; i++)
				;
			// search backward for a small entry
			for (j-- ; (j>=lo)&&(projection(j,axis,kdi)>pivot) ; j--)
				;
			perm_swap (kdi,i,j);
			}
	
		perm_swap (kdi,i,j);	// undo the last swap
		perm_swap (kdi,lo,j);	// move the pivot value to the proper location
	
		// warning: we must avoid returning j==hi (because build_kd_tree() would
		// recurse forever);  if j<hi, we had at least one value larger than the
		// pivot;  when j==hi the pivot was the max value;  if the range had
		// only two values we are assured that these two are sorted
	
		if      (j < hi)     return j;
		else if (hi-lo == 2) return hi-1;
	
		// otherwise, we need to partition them again, leaving out hi;  looping
		// back is equivalent to tail recursion, with this call:
		//    return partition_segments (lo, hi-1, axis, kdi);

		hi--;
		}

	}

//----------
//
// best_predecessor--
//	Find a segment's best predecessor chain.
//
// The best predecessor of a segment is the chain, among all those starting
// strictly before the segment in both x and y, that scores the highest when
// connected to this segment.
//
// This routine searches a subtree for such predecessor chains, pruning
// subtrees which contain no predecessor segments, or in which no segment can
// exceed the given lower bound.
//
//----------
//
// Arguments:
//	kdnode*		subtree:	The K-d (sub)tree to search.
//	int			axis:		Dimension/axis to partition on
//							  0 => diagonal (pos1 - pos2)
//							  1 => pos2
//	int			lowerBound:	Lower bound of chain score that must be achieved.
//	bestpred	bp:			The best predecessor found so far.
//	kdinfo*		kdi:		'Global' data.
//
// Returns:
//	The ending segment index (bestpred.num) and score (bestpred.contrib) of the
//	best predecessor chain.
//
//----------

#define debugChaining_4                                                      \
	if (chain_dbgChaining)                                                   \
		{                                                                    \
		if (bp.num == noPred)                                                \
			fprintf (stderr, "    bestwas=(none) scorewas=%.2f\n",           \
			                 bp.contrib);                                    \
		else                                                                 \
			fprintf (stderr, "    bestwas=%u scorewas=%.2f\n",               \
			                 bp.num, bp.contrib);                            \
		}

#define debugChaining_5                                                      \
	if (chain_dbgChaining)                                                   \
		{                                                                    \
		fprintf (stderr, "  cand=%u score=%.2f (from %.2f)\n",               \
		                 j, predScore, kdi->chainScore[j]);                  \
		}


static bestpred best_predecessor
   (kdnode*		subtree,
	int			axis,
	bigscore	lowerBound,
	bestpred	bp,
	const kdinfo* const kdi)
	{
	bigscore	predScore;
	u32			i, j;
	segment*	s;

	if (subtree == NULL)
		suicide ("(in best_predecessor, NULL subtree)");

	if (bp.contrib >= subtree->maxChainScore - lowerBound)
		return bp;

	if (!valid_kdnode(subtree))
		suicide ("(in best_predecessor, invalid subtree)");

	// if we're at a leaf, search over all segments in the leaf

	if (subtree->isBucket)
		{
		for (i=subtree->loIx ; i<=subtree->hiIx ; i++)
			{
			j = kdi->perm[i];   // kdi is the segment we want to add to the chain
			s = &kdi->seg[j];   // s is the candidate to be a predecessor
			if ((s->pos1 >= kdi->x) || (s->pos2 >= kdi->y))
				continue;
			predScore = kdi->chainScore[j] - kdi->connect(s, kdi->query, kdi->scale);
			debugChaining_4;
			if (predScore > bp.contrib) { bp.contrib = predScore;  bp.num = j; }
			debugChaining_5;
			}
		}

	// if we're at a node cut by y, search over both subtrees, pruning the high
	// subtree if all its segments have y greater than our query

	else if (axis == 1)
		{
		if (((sgnpos) kdi->y) >= subtree->cutVal)
			bp = best_predecessor (subtree->hiSon, lowerBound, 1-axis, bp, kdi);
		bp = best_predecessor (subtree->loSon, lowerBound, 1-axis, bp, kdi);
		}

	// if we're at a node cut by the diagonal, search over search both subtrees,
	// adjusting the lower bound accordingly
	// nota bene: diff>0 => query diagonal is below cut

	else // if (axis == 0)
		{
		bigscore diff = kdi->diag - subtree->cutVal;
		if (diff >= 0)	// query diagonal is southeast of (or same as) cut
			{
			bp = best_predecessor (subtree->hiSon, 1-axis, lowerBound,         bp, kdi);
			bp = best_predecessor (subtree->loSon, 1-axis, diff*kdi->diagPen,  bp, kdi);
			}
		else			// query diagonal is northwest of cut
			{
			bp = best_predecessor (subtree->loSon, 1-axis, lowerBound,         bp, kdi);
			bp = best_predecessor (subtree->hiSon, 1-axis, -diff*kdi->antiPen, bp, kdi);
			}
		}

	return bp;
	}

//----------
//
// propagate_max_score--
//	Propagate the best score for any chain ending at a particular segment to all
//	the (sub)trees that contain that segment.
//
//----------
//
// Arguments:
//	kdnode*		subtree:	The K-d (sub)tree to operate on.
//	bigscore	s:			The score to propagate.
//	u32			ix:			The index of the segment that has that score.
//
// Returns:
//	(nothing)
//
//----------
//
// Older, recursive version looked like this:
//
//	if (subtree != NULL)
//		{
//		if (s > subtree->maxChainScore) subtree->maxChainScore = s;
//		if (subtree->hiIx >= ix) propagate_max_score (subtree->loSon, s, ix);
//		                    else propagate_max_score (subtree->hiSon, s, ix);
//		}
//
//----------

static void propagate_max_score
   (kdnode*		subtree,
	bigscore	s,
	u32			ix)
	{
	while (subtree != NULL)
		{
		if (s > subtree->maxChainScore)
			subtree->maxChainScore = s;
		if (ix <= subtree->hiIx) subtree = subtree->loSon;
		                    else subtree = subtree->hiSon;
		}
	}

//----------
//
// chain_zero_stats--
//	Clear the statistics for this module.
//
//----------
//
// Arguments:
//	(none)
//
// Returns:
//	(nothing)
//
//----------

void chain_zero_stats
   (void)
	{
#ifdef collect_stats

	// set 'em en masse to zero

	memset (&chainStats, 0, sizeof(chainStats));

	// set any values that might be floating point to zero (fp bit pattern for
	// zero may not be all-bits-zero)

	// (none to set, yet)

#endif // collect_stats
	}

//----------
//
// chain_show_stats--
//	Show the statistics that have been collected for this module.
//
//----------
//
// Arguments:
//	FILE*		f:	The file to print the stats to.
//
// Returns:
//	(nothing)
//
//----------

void chain_show_stats
   (arg_dont_complain(FILE* f))
	{
#ifdef collect_stats
	if (f == NULL) return;
	fprintf (f, "# anchors to chain: %s\n", commatize(chainStats.numAnchors));
	fprintf (f, " segments in chain: %s\n", commatize(chainStats.numSegments));
	fprintf (f, "-------------------\n");
#endif // collect_stats
	}

void chain_generic_stats
   (arg_dont_complain(FILE* f),
    arg_dont_complain(void (*func) (FILE*, const char*, ...)))
	{
#ifdef collect_stats
	if (f == NULL) return;
	(*func) (f, "num_anchors=%d\n",   chainStats.numAnchors);
	(*func) (f, "num_segments=%d\n",  chainStats.numSegments);
#endif // collect_stats
	}

