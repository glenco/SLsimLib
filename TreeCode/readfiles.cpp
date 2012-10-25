

#include <slsimlib.h>

#define NR_END 1

PosType **PosTypeMatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a PosType matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	PosType **m;

	/* allocate pointers to rows */
	m=(PosType **) malloc((size_t)((nrow+NR_END)*sizeof(PosType*)));
	if (!m) {ERROR_MESSAGE(); std::cout << "ERROR: PosTypeMatrix" << std::endl << "allocation failure 1" << std::endl << std::endl; exit(0);}
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(PosType *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(PosType)));
	if (!m[nrl]) {ERROR_MESSAGE();  std::cout << "ERROR: PosTypeMatrix" << std::endl << "allocation failure 2" << std::endl << std::endl; exit(0);}
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}
void free_PosTypeMatrix(PosType **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((char*) (m[nrl]+ncl-NR_END));
	free((char*) (m+nrl-NR_END));
}

#undef NR_END
