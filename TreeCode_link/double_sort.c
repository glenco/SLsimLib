/* double_sort is adapted from NR sort2 to take an unsigned long */
/* second array brr */
/* note #undef's at end of file */
#include "Tree.h"
#define NRANSI
#include "../../Library/Recipes/nrutil.h"
//#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

void double_sort(unsigned long n, double *arr, unsigned long *brr)
{
  unsigned long i,ir=n,j,k,l=1,*istack,b;
  int jstack=0;
  double a,temp;
  
  istack=lvector(1,NSTACK);
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	a=arr[j];
	b=brr[j];
	for (i=j-1;i>=l;i--) {
	  if (arr[i] <= a) break;
	  arr[i+1]=arr[i];
	  brr[i+1]=brr[i];
	}
	arr[i+1]=a;
	brr[i+1]=b;
      }
      if (!jstack) {
	free_lvector(istack,1,NSTACK);
	return;
      }
      ir=istack[jstack];
      l=istack[jstack-1];
      jstack -= 2;
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1])
	SWAP(brr[k],brr[l+1])
	if (arr[l] > arr[ir]) {
	  SWAP(arr[l],arr[ir])
	    SWAP(brr[l],brr[ir])
	    }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	  SWAP(brr[l+1],brr[ir])
	  }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1])
	  SWAP(brr[l],brr[l+1])
	  }
      i=l+1;
      j=ir;
      a=arr[l+1];
      b=brr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j])
	  SWAP(brr[i],brr[j])
	  }
      arr[l+1]=arr[j];
      arr[j]=a;
      brr[l+1]=brr[j];
      brr[j]=b;
      jstack += 2;
      if (jstack > NSTACK) nrerror("NSTACK too small in double_sort");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
}

void double_sort_points(unsigned long n, double *arr, Point *brr){
  /* arr array uses NR standard indexing i.e arr[1...n] */
  /* but brr[0..n-1] */
  /* if the point array is two-way-coupled to another point array */
  /* the image pointers of that array will follow sort */
  /* if the array is not  two-way-coupled to another the image */
  /* pointers in the other array will be untouched */
  unsigned long i,ir=n,j,k,l=1,*istack;
  long jstack=0;
  double a,temp;
  Point b;

  istack=lvector(1,NSTACK);
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	a=arr[j];
	/*b=brr[j-1];*/
	 PointCopy(&b,&brr[j-1]);
	for (i=j-1;i>=l;i--) {
	  if (arr[i] <= a) break;
	  arr[i+1]=arr[i];
	  /*brr[i]=brr[i-1];*/
	   PointCopy(&brr[i],&brr[i-1]);
	}
	arr[i+1]=a;
	/*brr[i]=b;*/
	 PointCopy(&brr[i],&b);
      }
      if (!jstack) {
	free_lvector(istack,1,NSTACK);
	return;
      }
      ir=istack[jstack];

      l=istack[jstack-1];
      jstack -= 2;
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1])
	SwapPointsInArray(&brr[k-1],&brr[l]);
	if (arr[l] > arr[ir]) {
	  SWAP(arr[l],arr[ir])
	    SwapPointsInArray(&brr[l-1],&brr[ir-1]);
	    }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	  SwapPointsInArray(&brr[l],&brr[ir-1]);
	  }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1])
	  SwapPointsInArray(&brr[l-1],&brr[l]);
	  }
      i=l+1;
      j=ir;
      a=arr[l+1];
      /*b=brr[l];*/
       PointCopy(&b,&brr[l]);
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j])
	  SwapPointsInArray(&brr[i-1],&brr[j-1]);
	  }
      arr[l+1]=arr[j];
      arr[j]=a;
      /*brr[l]=brr[j-1];*/
       PointCopy(&brr[l],&brr[j-1]);
      /*brr[j-1]=b;*/
       PointCopy(&brr[j-1],&b);
      jstack += 2;
      if (jstack > NSTACK) nrerror("NSTACK too small in double_sort_points");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
}

#undef M
#undef NSTACK
#undef SWAP
#undef NRANSI


void quicksortPoints(Point *pointarray,double *arr,unsigned long N){
	double pivotvalue;
	unsigned long pivotindex,newpivotindex,i;
	void swap(double *a,double *b);

	if(N <= 1) return ;

	// pick pivot as the median of the first, last and middle values
	if ((arr[0] >= arr[N/2] && arr[0] <= arr[N-1])
			|| (arr[0] >= arr[N-1] && arr[0] <= arr[N/2])) pivotindex = 0;
	else if ((arr[N/2] >= arr[0] && arr[N/2] <= arr[N-1])
			|| (arr[N/2] >= arr[N-1] && arr[N/2] <= arr[0])) pivotindex = N/2;
	else pivotindex = N-1;
	pivotvalue=arr[pivotindex];

	// move pivet to end of array
	swap(&arr[pivotindex],&arr[N-1]);
	SwapPointsInArray(&pointarray[pivotindex],&pointarray[N-1]);
	newpivotindex=0;

	// partition list and array
	for(i=0;i<N;++i){
		if(arr[i] <= pivotvalue){
			swap(&arr[newpivotindex],&arr[i]);
			SwapPointsInArray(&pointarray[newpivotindex],&pointarray[i]);
			++newpivotindex;
		}
	}
	--newpivotindex;

	quicksortPoints(pointarray,arr,newpivotindex);
	quicksortPoints(&pointarray[newpivotindex+1],&arr[newpivotindex+1],N-newpivotindex-1);

	return ;
}

void quicksort(unsigned long *particles,double *arr,unsigned long N){
	double pivotvalue;
	unsigned long pivotindex,newpivotindex,i;
	void swap(double *a,double *b);
	void swapLong(unsigned long *a,unsigned long *b);

	if(N <= 1) return ;

	// pick pivot as the median of the first, last and middle values
	if ((arr[0] >= arr[N/2] && arr[0] <= arr[N-1])
			|| (arr[0] >= arr[N-1] && arr[0] <= arr[N/2])) pivotindex = 0;
	else if ((arr[N/2] >= arr[0] && arr[N/2] <= arr[N-1])
			|| (arr[N/2] >= arr[N-1] && arr[N/2] <= arr[0])) pivotindex = N/2;
	else pivotindex = N-1;
	pivotvalue=arr[pivotindex];

	// move pivet to end of array
	swap(&arr[pivotindex],&arr[N-1]);
	//SwapPointsInArray(&pointarray[pivotindex],&pointarray[N-1]);
	swapLong(&particles[pivotindex],&particles[N-1]);
	newpivotindex=0;

	// partition list and array
	for(i=0;i<N;++i){
		if(arr[i] <= pivotvalue){
			swap(&arr[newpivotindex],&arr[i]);
			//SwapPointsInArray(&pointarray[newpivotindex],&pointarray[i]);
			swapLong(&particles[newpivotindex],&particles[i]);
			++newpivotindex;
		}
	}
	--newpivotindex;

	quicksort(particles,arr,newpivotindex);
	quicksort(&particles[newpivotindex+1],&arr[newpivotindex+1],N-newpivotindex-1);

	return ;
}

/*
 * Partitions arr[] and particles[] into those with x <= pivotvalue and those with
 *  x > pivotvalue  pivotindex is left at first array value with x > pivotvalue
 */
void quickPartition(double pivotvalue,unsigned long *pivotindex,unsigned long *particles
		,double *arr,unsigned long N){
	unsigned long i;
	void swap(double *a,double *b);
	void swapLong(unsigned long *a,unsigned long *b);

	*pivotindex=0;

	for(i=0;i<N;++i){
		if(arr[i] <= pivotvalue){
			swap(&arr[*pivotindex],&arr[i]);
			swapLong(&particles[*pivotindex],&particles[i]);
			++(*pivotindex);
		}
	}

	return ;
}
void quickPartitionPoints(double pivotvalue,unsigned long *pivotindex
		,Point *pointarray,double *arr,unsigned long N){
	unsigned long i;
	void swap(double *a,double *b);
	void swapLong(unsigned long *a,unsigned long *b);

	*pivotindex=0;

	for(i=0;i<N;++i){
		if(arr[i] <= pivotvalue){
			swap(&arr[*pivotindex],&arr[i]);
			SwapPointsInArray(&pointarray[*pivotindex],&pointarray[i]);
			++(*pivotindex);
		}
	}

	return ;
}

void swap(double *a,double *b){
	double tmp;
	tmp=*a;
	*a=*b;
	*b=tmp;
}
void swapLong(unsigned long *a,unsigned long *b){
	unsigned long tmp;
	tmp=*a;
	*a=*b;
	*b=tmp;
}
