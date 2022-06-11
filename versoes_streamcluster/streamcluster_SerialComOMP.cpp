#include <assert.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>

#include <fstream>
#include <iostream>

#ifdef ENABLE_PARSEC_HOOKS
#include <hooks.h>
#endif

using namespace std;

#define MAXNAMESIZE 1024  // max filename length
#define SEED 1
#define SP 1	// number of repetitions of speedy must be >=1
#define ITER 3	// iterate ITER* k log k times; ITER >= 1

#define CACHE_LINE 32  // cache line in byte

typedef struct {
	float weight;
	float *coord;
	long assign; /* number of point where this one is assigned */
	float cost;	 /* cost of that assignment, weight*distance */
} Point;

typedef struct {
	long num; /* number of points; may not be N if this is a sample */
	int dim;  /* dimensionality */
	Point *p; /* the array itself */
} Points;

static bool *switch_membership;	 //whether to switch membership in pgain
static bool *is_center;			 //whether a point is a center
static int *center_table;		 //index table of centers

static int nproc;  //# of threads
//int gg = 0;
//double* timesgg;
//double ggsum = 0;

float dist(Point p1, Point p2, int dim);

int isIdentical(float *i, float *j, int D) {
	int a = 0;
	int equal = 1;

	while (equal && a < D) {
		if (i[a] != j[a])
			equal = 0;
		else
			a++;
	}
	if (equal)
		return 1;
	else
		return 0;
}

static int floatcomp(const void *i, const void *j) {
	float a, b;
	a = *(float *)(i);
	b = *(float *)(j);
	if (a > b) return (1);
	if (a < b) return (-1);
	return (0);
}

void shuffle(Points *points) {
	long i, j;
	Point temp;
	for (i = 0; i < points->num - 1; i++) {
		j = (lrand48() % (points->num - i)) + i;
		temp = points->p[i];
		points->p[i] = points->p[j];
		points->p[j] = temp;
	}
}

void intshuffle(int *intarray, int length) {
	long i, j;
	int temp;
	for (i = 0; i < length; i++) {
		j = (lrand48() % (length - i)) + i;
		temp = intarray[i];
		intarray[i] = intarray[j];
		intarray[j] = temp;
	}
}

float dist(Point p1, Point p2, int dim) {
	int i;
	float result = 0.0;
	for (i = 0; i < dim; i++)
		result += (p1.coord[i] - p2.coord[i]) * (p1.coord[i] - p2.coord[i]);
	return (result);
}

float pspeedy(Points *points, float z, long *kcenter) {
	long k1 = 0;
	long k2 = points->num;

	static double totalcost;

	static bool open = false;
	static double costs;  //cost for each thread.
	static int i;

	for (int k = k1; k < k2; k++) {
		float distance = dist(points->p[k], points->p[0], points->dim);
		points->p[k].cost = distance * points->p[k].weight;
		points->p[k].assign = 0;
	}

	*kcenter = 1;

	for (i = 1; i < points->num; i++) {
		bool to_open = ((float)lrand48() / (float)INT_MAX) < (points->p[i].cost / z);
		if (to_open) {
			(*kcenter)++;

			open = true;

			for (int k = k1; k < k2; k++) {
				float distance = dist(points->p[i], points->p[k], points->dim);
				if (distance * points->p[k].weight < points->p[k].cost) {
					points->p[k].cost = distance * points->p[k].weight;
					points->p[k].assign = i;
				}
			}

			open = false;
		}
	}

	open = false;
	double mytotal = 0;
	for (int k = k1; k < k2; k++) {
		mytotal += points->p[k].cost;
	}
	costs = mytotal;

	// aggregate costs from each thread
	totalcost = z * (*kcenter);
	totalcost += costs;

	return (totalcost);
}

double pgain(long x, Points *points, double z, long int *numcenters, int pid)
{
	// my block
	long bsize = points->num / nproc;
	long k1 = bsize * pid;
	long k2 = k1 + bsize;
	if (pid == nproc - 1)
		k2 = points->num;

	int i;
	int number_of_centers_to_close = 0;

	static double *work_mem;
	static double gl_cost_of_opening_x;
	static int gl_number_of_centers_to_close;

	// each thread takes a block of working_mem.
	int stride = *numcenters + 2;
	// make stride a multiple of CACHE_LINE
	int cl = CACHE_LINE / sizeof(double);
	if (stride % cl != 0)
	{
		stride = cl * (stride / cl + 1);
	}
	int K = stride - 2; // K==*numcenters

	// my own cost of opening x
	double cost_of_opening_x = 0;

	if (pid == 0)
	{
		work_mem = (double *)malloc(stride * (nproc + 1) * sizeof(double));
		gl_cost_of_opening_x = 0;
		gl_number_of_centers_to_close = 0;
	}

#pragma omp barrier

	/*For each center, we have a *lower* field that indicates
	how much we will save by closing the center.
	Each thread has its own copy of the *lower* fields as an array.
	We first build a table to index the positions of the *lower* fields.
  */
#pragma omp parallel
	{
		int count = 0;
		int i;
#pragma omp for private(i)
		for (i = k1; i < k2; i++)
		{
			if (is_center[i])
			{
#pragma omp critical
				center_table[i] = count++;
			}
		}
		work_mem[pid * stride] = count;

#pragma omp barrier

		if (pid == 0)
		{
			int accum = 0;
			int p;
#pragma omp for private(p, stride)
			for (p = 0; p < nproc; p++)
			{
				int tmp = (int)work_mem[p * stride];
				work_mem[p * stride] = accum;
#pragma omp atomic
				accum += tmp;
			}
		}

#pragma omp barrier

#pragma omp for private(i, pid, stride)
		for (int i = k1; i < k2; i++)
		{
			if (is_center[i])
			{
#pragma omp critical
				center_table[i] += (int)work_mem[pid * stride];
			}
		}

// now we finish building the table. clear the working memory.
#pragma omp sections
		{
#pragma omp section
			memset(switch_membership + k1, 0, (k2 - k1) * sizeof(bool));

#pragma omp section
			memset(work_mem + pid * stride, 0, stride * sizeof(double));
		}
		if (pid == 0)
			memset(work_mem + nproc * stride, 0, stride * sizeof(double));

#pragma omp barrier

		// my *lower* fields
		double *lower = &work_mem[pid * stride];
		// global *lower* fields
		double *gl_lower = &work_mem[nproc * stride];

		float x_cost;
		float current_cost;
#pragma omp for private(i, x_cost, current_cost)
		for (i = k1; i < k2; i++)
		{
			x_cost = dist(points->p[i], points->p[x], points->dim) * points->p[i].weight;
			current_cost = points->p[i].cost;

			if (x_cost < current_cost)
			{
				// point i would save cost just by switching to x
				// (note that i cannot be a median,
				// or else dist(p[i], p[x]) would be 0)

				switch_membership[i] = 1;
#pragma omp critical
				cost_of_opening_x += x_cost - current_cost;
			}
			else
			{
				// cost of assigning i to x is at least current assignment cost of i

				// consider the savings that i's **current** median would realize
				// if we reassigned that median and all its members to x;
				// note we've already accounted for the fact that the median
				// would save z by closing; now we have to subtract from the savings
				// the extra cost of reassigning that median and its members
				int assign = points->p[i].assign;
#pragma omp critical
				lower[center_table[assign]] += current_cost - x_cost;
			}
		}

#pragma omp barrier

		// at this time, we can calculate the cost of opening a center
		// at x; if it is negative, we'll go through with opening it
		double low;
#pragma omp for private(i, low)
		for (int i = k1; i < k2; i++)
		{
			if (is_center[i])
			{
				low = z;
				// aggregate from all threads
				int p;
				for (p = 0; p < nproc; p++)
				{
					low += work_mem[center_table[i] + p * stride];
				}
#pragma omp critical
				gl_lower[center_table[i]] = low;
				if (low > 0)
				{
// i is a median, and
// if we were to open x (which we still may not) we'd close i

// note, we'll ignore the following quantity unless we do open x
#pragma omp atomic
					++number_of_centers_to_close;
#pragma omp critical
					cost_of_opening_x -= low;
				}
			}
		}

		work_mem[pid * stride + K] = number_of_centers_to_close;
		work_mem[pid * stride + K + 1] = cost_of_opening_x;

#pragma omp barrier
		//  printf("thread %d cost complete\n",pid);
		if (pid == 0)
		{
			gl_cost_of_opening_x = z;
			// aggregate
			int p;
#pragma omp for private(p)
			for (p = 0; p < nproc; p++)
			{
				gl_number_of_centers_to_close += (int)work_mem[p * stride + K];
				gl_cost_of_opening_x += work_mem[p * stride + K + 1];
			}
		}

#pragma omp barrier
		// Now, check whether opening x would save cost; if so, do it, and
		// otherwise do nothing

		if (gl_cost_of_opening_x < 0)
		{
			//  we'd save money by opening x; we'll do it
			bool close_center;
#pragma omp for private(i, close_center)
			for (int i = k1; i < k2; i++)
			{
				close_center = gl_lower[center_table[points->p[i].assign]] > 0;
				if (switch_membership[i] || close_center)
				{
					// Either i's median (which may be i itself) is closing,
					// or i is closer to x than to its current median
					points->p[i].cost = points->p[i].weight *
										dist(points->p[i], points->p[x], points->dim);
					points->p[i].assign = x;
				}
			}
#pragma omp for private(i)
			for (i = k1; i < k2; i++)
			{
				if (is_center[i] && gl_lower[center_table[i]] > 0)
				{
					is_center[i] = false;
				}
			}
			if (x >= k1 && x < k2)
			{
				is_center[x] = true;
			}

			if (pid == 0)
			{
				*numcenters = *numcenters + 1 - gl_number_of_centers_to_close;
			}
		}
		else
		{
			if (pid == 0)
				gl_cost_of_opening_x = 0; // the value we'll return
		}
#pragma omp barrier
	} // parallel
	if (pid == 0)
	{
		free(work_mem);
		//    free(is_center);
		//    free(switch_membership);
		//    free(proc_cost_of_opening_x);
		//    free(proc_number_of_centers_to_close);
	}

	return -gl_cost_of_opening_x;
}

float pFL(Points *points, int *feasible, int numfeasible,
		  float z, long *k, double cost, long iter, float e) {
	long i;
	long x;
	double change;
	long numberOfPoints;

	change = cost;

	while (change / cost > 1.0 * e) {
		change = 0.0;
		numberOfPoints = points->num;
		/* randomize order in which centers are considered */
		intshuffle(feasible, numfeasible);

		for (i = 0; i < iter; i++) {
			x = i % numfeasible;
			change += pgain(feasible[x], points, z, k);
		}
		cost -= change;
	}
	return (cost);
}

int selectfeasible_fast(Points *points, int **feasible, int kmin) {
	int numfeasible = points->num;
	if (numfeasible > (ITER * kmin * log((double)kmin)))
		numfeasible = (int)(ITER * kmin * log((double)kmin));
	*feasible = (int *)malloc(numfeasible * sizeof(int));

	float *accumweight;
	float totalweight;

	long k1 = 0;
	long k2 = numfeasible;

	float w;
	int l, r, k;

	if (numfeasible == points->num) {
		for (int i = k1; i < k2; i++)
			(*feasible)[i] = i;
		return numfeasible;
	}

	accumweight = (float *)malloc(sizeof(float) * points->num);

	accumweight[0] = points->p[0].weight;
	totalweight = 0;
	for (int i = 1; i < points->num; i++) {
		accumweight[i] = accumweight[i - 1] + points->p[i].weight;
	}
	totalweight = accumweight[points->num - 1];

	for (int i = k1; i < k2; i++) {
		w = (lrand48() / (float)INT_MAX) * totalweight;
		l = 0;
		r = points->num - 1;
		if (accumweight[0] > w) {
			(*feasible)[i] = 0;
			continue;
		}
		while (l + 1 < r) {
			k = (l + r) / 2;
			if (accumweight[k] > w) {
				r = k;
			} else {
				l = k;
			}
		}
		(*feasible)[i] = r;
	}

	free(accumweight);

	return numfeasible;
}

float pkmedian(Points *points, long kmin, long kmax, long *kfinal) {
	int i;
	double cost;
	double lastcost;
	double hiz, loz, z;

	static long k;
	static int *feasible;
	static int numfeasible;
	static double hizs;

	hizs = 0.0;
	hiz = loz = 0.0;
	long numberOfPoints = points->num;
	long ptDimension = points->dim;

	//my block
	long k1 = 0;
	long k2 = points->num;

	double myhiz = 0;
	for (long kk = k1; kk < k2; kk++) {
		myhiz += dist(points->p[kk], points->p[0],
					  ptDimension) *
				 points->p[kk].weight;
	}
	hizs = myhiz;
	hiz += hizs;

	loz = 0.0;
	z = (hiz + loz) / 2.0;
	/* NEW: Check whether more centers than points! */
	if (points->num <= kmax) {
		/* just return all points as facilities */
		for (long kk = k1; kk < k2; kk++) {
			points->p[kk].assign = kk;
			points->p[kk].cost = 0;
		}
		cost = 0;
		*kfinal = k;

		return cost;
	}

	shuffle(points);
	cost = pspeedy(points, z, &k);

	i = 0;
	/* give speedy SP chances to get at least kmin/2 facilities */
	while ((k < kmin) && (i < SP)) {
		cost = pspeedy(points, z, &k);
		i++;
	}

	/* if still not enough facilities, assume z is too high */
	while (k < kmin) {
		if (i >= SP) {
			hiz = z;
			z = (hiz + loz) / 2.0;
			i = 0;
		}
		shuffle(points);
		cost = pspeedy(points, z, &k);
		i++;
	}

	numfeasible = selectfeasible_fast(points, &feasible, kmin);
	for (int i = 0; i < points->num; i++) {
		is_center[points->p[i].assign] = true;
	}

	while (1) {
		/* first get a rough estimate on the FL solution */
		lastcost = cost;
		cost = pFL(points, feasible, numfeasible,
				   z, &k, cost, (long)(ITER * kmax * log((double)kmax)), 0.1);

		/* if number of centers seems good, try a more accurate FL */
		if (((k <= (1.1) * kmax) && (k >= (0.9) * kmin)) ||
			((k <= kmax + 2) && (k >= kmin - 2))) {
			cost = pFL(points, feasible, numfeasible,
					   z, &k, cost, (long)(ITER * kmax * log((double)kmax)), 0.001);
		}

		if (k > kmax) {
			loz = z;
			z = (hiz + loz) / 2.0;
			cost += (z - loz) * k;
		}
		if (k < kmin) {
			hiz = z;
			z = (hiz + loz) / 2.0;
			cost += (z - hiz) * k;
		}

		if (((k <= kmax) && (k >= kmin)) || ((loz >= (0.999) * hiz))) {
			break;
		}
	}

	*kfinal = k;

	return cost;
}

int contcenters(Points *points) {
	long i, ii;
	float relweight;

	for (i = 0; i < points->num; i++) {
		/* compute relative weight of this point to the cluster */
		if (points->p[i].assign != i) {
			relweight = points->p[points->p[i].assign].weight + points->p[i].weight;
			relweight = points->p[i].weight / relweight;
			for (ii = 0; ii < points->dim; ii++) {
				points->p[points->p[i].assign].coord[ii] *= 1.0 - relweight;
				points->p[points->p[i].assign].coord[ii] +=
					points->p[i].coord[ii] * relweight;
			}
			points->p[points->p[i].assign].weight += points->p[i].weight;
		}
	}

	return 0;
}

void copycenters(Points *points, Points *centers, long *centerIDs, long offset) {
	long i;
	long k;

	bool *is_a_median = (bool *)calloc(points->num, sizeof(bool));

	/* mark the centers */
	for (i = 0; i < points->num; i++) {
		is_a_median[points->p[i].assign] = 1;
	}

	k = centers->num;

	/* count how many  */
	for (i = 0; i < points->num; i++) {
		if (is_a_median[i]) {
			memcpy(centers->p[k].coord, points->p[i].coord, points->dim * sizeof(float));
			centers->p[k].weight = points->p[i].weight;
			centerIDs[k] = i + offset;
			k++;
		}
	}

	centers->num = k;

	free(is_a_median);
}

struct pkmedian_arg_t {
	Points *points;
	long kmin;
	long kmax;
	long *kfinal;
};

void *localSearchSub(void *arg_) {
	pkmedian_arg_t *arg = (pkmedian_arg_t *)arg_;
	pkmedian(arg->points, arg->kmin, arg->kmax, arg->kfinal);

	return NULL;
}

void localSearch(Points *points, long kmin, long kmax, long *kfinal) {
	pkmedian_arg_t arg;

	arg.points = points;
	arg.kmin = kmin;
	arg.kmax = kmax;
	arg.kfinal = kfinal;

	localSearchSub(&arg);
}

class PStream {
   public:
	virtual size_t read(float *dest, int dim, int num) = 0;
	virtual int ferror() = 0;
	virtual int feof() = 0;
	virtual ~PStream() {
	}
};

class SimStream : public PStream {
   public:
	SimStream(long n_) {
		n = n_;
	}
	size_t read(float *dest, int dim, int num) {
		size_t count = 0;
		for (int i = 0; i < num && n > 0; i++) {
			for (int k = 0; k < dim; k++) {
				dest[i * dim + k] = lrand48() / (float)INT_MAX;
			}
			n--;
			count++;
		}
		return count;
	}
	int ferror() {
		return 0;
	}
	int feof() {
		return n <= 0;
	}
	~SimStream() {
	}

   private:
	long n;
};

class FileStream : public PStream {
   public:
	FileStream(char *filename) {
		fp = fopen(filename, "rb");
		if (fp == NULL) {
			fprintf(stderr, "error opening file %s\n.", filename);
			exit(1);
		}
	}
	size_t read(float *dest, int dim, int num) {
		return std::fread(dest, sizeof(float) * dim, num, fp);
	}
	int ferror() {
		return std::ferror(fp);
	}
	int feof() {
		return std::feof(fp);
	}
	~FileStream() {
		fprintf(stderr, "closing file stream\n");
		fclose(fp);
	}

   private:
	FILE *fp;
};

void outcenterIDs(Points *centers, long *centerIDs, char *outfile) {
	FILE *fp = fopen(outfile, "w");
	if (fp == NULL) {
		fprintf(stderr, "error opening %s\n", outfile);
		exit(1);
	}
	int *is_a_median = (int *)calloc(sizeof(int), centers->num);
	for (int i = 0; i < centers->num; i++) {
		is_a_median[centers->p[i].assign] = 1;
	}

	for (int i = 0; i < centers->num; i++) {
		if (is_a_median[i]) {
			fprintf(fp, "%u\n", centerIDs[i]);
			fprintf(fp, "%lf\n", centers->p[i].weight);
			for (int k = 0; k < centers->dim; k++) {
				fprintf(fp, "%lf ", centers->p[i].coord[k]);
			}
			fprintf(fp, "\n\n");
		}
	}
	fclose(fp);
}

void streamCluster(PStream *stream,
				   long kmin, long kmax, int dim,
				   long chunksize, long centersize, char *outfile) {
	float *block = (float *)malloc(chunksize * dim * sizeof(float));
	float *centerBlock = (float *)malloc(centersize * dim * sizeof(float));
	long *centerIDs = (long *)malloc(centersize * dim * sizeof(long));

	if (block == NULL) {
		fprintf(stderr, "not enough memory for a chunk!\n");
		exit(1);
	}

	Points points;
	points.dim = dim;
	points.num = chunksize;
	points.p = (Point *)malloc(chunksize * sizeof(Point));

	for (int i = 0; i < chunksize; i++) {
		points.p[i].coord = &block[i * dim];
	}

	Points centers;
	centers.dim = dim;
	centers.p =

		(Point *)malloc(centersize * sizeof(Point));
	centers.num = 0;

	for (int i = 0; i < centersize; i++) {
		centers.p[i].coord = &centerBlock[i * dim];
		centers.p[i].weight = 1.0;
	}

	long IDoffset = 0;
	long kfinal;
	while (1) {
		size_t numRead = stream->read(block, dim, chunksize);
		//fprintf(stderr,"read %d points\n",numRead);

		if (stream->ferror() || numRead < (unsigned int)chunksize && !stream->feof()) {
			fprintf(stderr, "error reading data!\n");
			exit(1);
		}

		points.num = numRead;
		for (int i = 0; i < points.num; i++) {
			points.p[i].weight = 1.0;
		}

		switch_membership = (bool *)malloc(points.num * sizeof(bool));
		is_center = (bool *)calloc(points.num, sizeof(bool));
		center_table = (int *)malloc(points.num * sizeof(int));

		localSearch(&points, kmin, kmax, &kfinal);	// parallel
		contcenters(&points);						/* sequential */

		if (kfinal + centers.num > centersize) {
			//here we don't handle the situation where # of centers gets too large.
			fprintf(stderr, "oops! no more space for centers\n");
			exit(1);
		}

		copycenters(&points, &centers, centerIDs, IDoffset); /* sequential */
		IDoffset += numRead;

		free(is_center);
		free(switch_membership);
		free(center_table);

		if (stream->feof()) {
			break;
		}
	}

	//finally cluster all temp centers

	switch_membership = (bool *)malloc(centers.num * sizeof(bool));
	is_center = (bool *)calloc(centers.num, sizeof(bool));
	center_table = (int *)malloc(centers.num * sizeof(int));
	//timesgg = (double*)malloc(2000*sizeof(double));

	localSearch(&centers, kmin, kmax, &kfinal);	 // parallel
	contcenters(&centers);
	outcenterIDs(&centers, centerIDs, outfile);
}

int main(int argc, char **argv) {
	char *outfilename = new char[MAXNAMESIZE];
	char *infilename = new char[MAXNAMESIZE];
	long kmin, kmax, n, chunksize, clustersize;
	int dim;

#ifdef PARSEC_VERSION
#define __PARSEC_STRING(x) #x
#define __PARSEC_XSTRING(x) __PARSEC_STRING(x)
	fprintf(stderr, "PARSEC Benchmark Suite Version "__PARSEC_XSTRING(PARSEC_VERSION) "\n");
	fflush(NULL);
#else
	fprintf(stderr, "PARSEC Benchmark Suite\n");
	fflush(NULL);
#endif	//PARSEC_VERSION
#ifdef ENABLE_PARSEC_HOOKS
	__parsec_bench_begin(__parsec_streamcluster);
#endif

	if (argc < 10) {
		fprintf(stderr, "usage: %s k1 k2 d n chunksize clustersize infile outfile nproc\n",
				argv[0]);
		fprintf(stderr, "  k1:          Min. number of centers allowed\n");
		fprintf(stderr, "  k2:          Max. number of centers allowed\n");
		fprintf(stderr, "  d:           Dimension of each data point\n");
		fprintf(stderr, "  n:           Number of data points\n");
		fprintf(stderr, "  chunksize:   Number of data points to handle per step\n");
		fprintf(stderr, "  clustersize: Maximum number of intermediate centers\n");
		fprintf(stderr, "  infile:      Input file (if n<=0)\n");
		fprintf(stderr, "  outfile:     Output file\n");
		fprintf(stderr, "  nproc:       Number of threads to use\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "if n > 0, points will be randomly generated instead of reading from infile.\n");
		exit(1);
	}

	kmin = atoi(argv[1]);
	kmax = atoi(argv[2]);
	dim = atoi(argv[3]);
	n = atoi(argv[4]);
	chunksize = atoi(argv[5]);
	clustersize = atoi(argv[6]);
	strcpy(infilename, argv[7]);
	strcpy(outfilename, argv[8]);
	nproc = atoi(argv[9]);
    omp_set_num_threads(nproc);
	srand48(SEED);
	PStream *stream;
	if (n > 0) {
		stream = new SimStream(n);
	} else {
		stream = new FileStream(infilename);
	}

#ifdef ENABLE_PARSEC_HOOKS
	__parsec_roi_begin();
#endif

	streamCluster(stream, kmin, kmax, dim, chunksize, clustersize, outfilename);

#ifdef ENABLE_PARSEC_HOOKS
	__parsec_roi_end();
#endif

	delete stream;
	
#ifdef ENABLE_PARSEC_HOOKS
	__parsec_bench_end();
#endif

	return 0;
}