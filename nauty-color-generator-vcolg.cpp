// nauty-color-generator-vcolg.cpp
// adapted from
/* vcolg.c version 2.0; B D McKay, May 11, 2017 */

#include "nauty-color-generator-vcolg.h"

#include "gtools.h"
#include "naugroup.h"
#include "nautinv.h"

#include <array>

static const TLS_ATTR CoffeeCode::NautyLink::VColGCallbackType *outproc;

static TLS_ATTR nauty_counter vc_nin, vc_nout;
static TLS_ATTR FILE *outfile;

#define MAXNV MAXN

static TLS_ATTR int col[MAXNV];
static TLS_ATTR boolean first;
static TLS_ATTR int lastreject[MAXNV];
static TLS_ATTR boolean lastrejok;
static TLS_ATTR unsigned long groupsize;
static TLS_ATTR unsigned long newgroupsize;
static TLS_ATTR boolean Tswitch;

static TLS_ATTR int fail_level;


/**************************************************************************/

static int
ismax(int *p, int n)
/* test if col^p <= col */
{
	int i, k;
	int fail;

	fail = 0;
	for (i = 0; i < n; ++i)
	{
		k = p[i];
		if (k > fail) fail = k;
		if (col[k] > col[i])
		{
			fail_level = fail;
			return FALSE;
		}
		else if (col[k] < col[i]) return TRUE;
	}

	++newgroupsize;
	return TRUE;
}

/**************************************************************************/

static void
testmax(int *p, int n, int *abort)
/* Called by allgroup2. */
{
	int i;

	if (first)
	{                       /* only the identity */
		first = FALSE;
		return;
	}

	if (!ismax(p, n))
	{
		*abort = 1;
		for (i = 0; i < n; ++i) lastreject[i] = p[i];
		lastrejok = TRUE;
	}
}

/**************************************************************************/

static int
trythisone(grouprec *group, graph *g, boolean digraph, int m, int n)
/* Try one solution, accept if maximal. */
/* Return value is level to return to. */
{
	int i, j;
	boolean accept;
	graph *gi;
	size_t ne;

	newgroupsize = 1;

	if (!group || groupsize == 1)
		accept = TRUE;
	else if (lastrejok && !ismax(lastreject, n))
		accept = FALSE;
	else if (lastrejok && groupsize == 2)
		accept = TRUE;
	else
	{
		newgroupsize = 1;
		first = TRUE;

		if (allgroup2(group, testmax) == 0)
			accept = TRUE;
		else
			accept = FALSE;
	}

	if (accept)
	{
		++vc_nout;

		(*outproc)(col, vc_nout);

		return n - 1;
	}
	else
		return fail_level - 1;
}

/**************************************************************************/

static int
scan(int level, graph *g, boolean digraph, int *prev, long minedges, long maxedges,
	long sofar, long numcols, grouprec *group, int m, int n)
	/* Recursive scan for default case */
	/* Returned value is level to return to. */
{
	int left;
	long min, max, k, ret;

	if (level == K_SYS)
		return trythisone(group, g, digraph, m, K_SYS);

	left = n - level - 1;
	min = minedges - sofar - numcols * left;
	if (min < 0) min = 0;
	max = maxedges - sofar;
	if (max >= numcols) max = numcols - 1;
	if (prev[level] >= 0 && col[prev[level]] < max)
		max = col[prev[level]];

	for (k = min; k <= max; ++k)
	{
		col[level] = k;
		ret = scan(level + 1, g, digraph, prev, minedges, maxedges, sofar + k, numcols, group, m, n);
		if (ret < level) return ret;
	}

	return level - 1;
}

/**************************************************************************/

static void
colourgraph(graph *g, int nfixed, long minedges, long maxedges,
	long numcols, int m, int n)
{
	static DEFAULTOPTIONS_GRAPH(options);
	statsblk stats;
	setword workspace[MAXNV];
	grouprec *group;
	int i, j, k, nloops;
	set *gi, *gj;
	int lab[MAXNV], ptn[MAXNV], orbits[MAXNV];
	boolean loop[MAXNV];
	int prev[MAXNV]; /* If >= 0, earlier point that must have greater colour */
	int weight[MAXNV];
	int region, start, stop;

	//if (n > MAXNV) gt_abort(">E vcolg: MAXNV exceeded\n");
	nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);

	nloops = 0;
	for (i = 0, gi = g; i < n; ++i, gi += m)
		if (ISELEMENT(gi, i))
		{
			DELELEMENT(gi, i);
			loop[i] = TRUE;
			++nloops;
		}
		else
			loop[i] = FALSE;

	for (region = 0; region < 2; ++region)
	{
		if (region == 0)
		{
			if (nfixed == 0) continue;
			start = 0;
			stop = nfixed;
			if (stop > n) stop = n;
		}
		else
		{
			if (nfixed >= n) continue;
			start = nfixed;
			stop = n;
		}

		for (i = start, gi = g + m * (size_t)start; i < stop; ++i, gi += m)
		{
			/* Find most recent equivalent j. */
			for (j = i - 1, gj = gi - m; j >= start; --j, gj -= m)
			{
				if (loop[j] != loop[i]) continue;
				for (k = 0; k < m; ++k) if (gi[k] != gj[k]) break;
				if (k < m)
				{
					FLIPELEMENT(gi, i); FLIPELEMENT(gj, j);
					for (k = 0; k < m; ++k) if (gi[k] != gj[k]) break;
					FLIPELEMENT(gi, i); FLIPELEMENT(gj, j);
				}
				if (k == m) break;
			}
			if (j >= start)
			{
				prev[i] = j;
				weight[i] = weight[j] + 1;
			}
			else
			{
				prev[i] = -1;
				weight[i] = 0;
			}
		}
	}

	if (n == 0)
	{
		scan(0, g, FALSE, prev, minedges, maxedges, 0, numcols, FALSE, m, n);
		return;
	}

	for (i = nfixed; i < n; ++i) weight[i] += nfixed;

	if (maxedges == NOLIMIT || maxedges > n*numcols) maxedges = n * numcols;
	if (n*numcols < minedges) return;

	options.userautomproc = groupautomproc;
	options.userlevelproc = grouplevelproc;
	options.defaultptn = FALSE;
	options.digraph = (nloops > 0);

	setlabptn(weight, lab, ptn, n);

	if (nloops > 0)
		for (i = 0, gi = g; i < n; ++i, gi += m)
			if (loop[i]) ADDELEMENT(gi, i);

	nauty(g, lab, ptn, NULL, orbits, &options, &stats, workspace, MAXNV, m, n, NULL);

	if (stats.grpsize2 == 0)
		groupsize = stats.grpsize1 + 0.1;
	else
		groupsize = 0;

	group = groupptr(FALSE);
	makecosetreps(group);

	if (stats.numorbits < n)
	{
		j = n;
		for (i = 0; i < n; ++i)
			if (orbits[i] < i && orbits[i] < j) j = orbits[i];

		for (i = j + 1; i < n; ++i)
			if (orbits[i] == j) prev[i] = j;
	}

	lastrejok = FALSE;
	for (i = 0; i < n; ++i) col[i] = 0;

	scan(0, g, FALSE, prev, minedges, maxedges, 0, numcols, group, m, n);
}

namespace CoffeeCode::NautyLink {

	// cast graph pointer type
	void vcolg(const VColGCallbackType& callback, void *g, size_t numcols)
	{
		outproc = &callback;
		colourgraph(reinterpret_cast<graph*>(g), K_SYS, 0, NOLIMIT, numcols, MAXM, MAXN);
	}
}