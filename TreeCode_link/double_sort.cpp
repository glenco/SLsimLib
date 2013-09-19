
/* double_sort is adapted from NR sort2 to take an unsigned long */
/* second array brr */
/* note #undef's at end of file */
#include "slsimlib.h"

#define NRANSI
#define M 7
#define NSTACK 50
namespace Utilities{
void double_sort(unsigned long n, double *arr, unsigned long *brr)
{
  unsigned long i,ir=n,j,k,l=1,*istack,b;
  int jstack=0;
  double a;
  
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
      std::swap(arr[k],arr[l+1]);
      std::swap(brr[k],brr[l+1]);
	if (arr[l] > arr[ir]) {
	  std::swap(arr[l],arr[ir]);
	  std::swap(brr[l],brr[ir]);
	    }
      if (arr[l+1] > arr[ir]) {
	std::swap(arr[l+1],arr[ir]);
	std::swap(brr[l+1],brr[ir]);
	  }
      if (arr[l] > arr[l+1]) {
	std::swap(arr[l],arr[l+1]);
	std::swap(brr[l],brr[l+1]);
	  }
      i=l+1;
      j=ir;
      a=arr[l+1];
      b=brr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	std::swap(arr[i],arr[j]);
	std::swap(brr[i],brr[j]);
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

/** \brief Sorts points in a point array.
 *
 *  arr array uses NR standard indexing i.e arr[1...n]
* but brr[0..n-1]
* if the point array is two-way-coupled to another point array
* the image pointers of that array will follow sort
* if the array is not  two-way-coupled to another the image
* pointers in the other array will be untouched
*/
void double_sort_points(unsigned long n, double *arr, Point *brr){
  unsigned long i,ir=n,j,k,l=1,*istack;
  long jstack=0;
  double a;
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
      std::swap(arr[k],arr[l+1]);
      SwapPointsInArray(&brr[k-1],&brr[l]);
      if (arr[l] > arr[ir]) {
    	  std::swap(arr[l],arr[ir]);
	      SwapPointsInArray(&brr[l-1],&brr[ir-1]);
      }
      if (arr[l+1] > arr[ir]) {
    	  std::swap(arr[l+1],arr[ir]);
    	  SwapPointsInArray(&brr[l],&brr[ir-1]);
	  }
      if (arr[l] > arr[l+1]) {
    	  std::swap(arr[l],arr[l+1]);
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
    	  std::swap(arr[i],arr[j]);
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
#undef NRANSI

/** \ingroup Utill
 * \brief Sorts points from smallest to largest according to the value of arr[].
 * Sorts arr[] and pointarray[] simultaneously.
 */

void quicksortPoints(Point *pointarray,double *arr,unsigned long N){
	double pivotvalue;
	unsigned long pivotindex,newpivotindex,i;
	
	if(N <= 1) return ;

	// pick pivot as the median of the first, last and middle values
	if ((arr[0] >= arr[N/2] && arr[0] <= arr[N-1])
			|| (arr[0] >= arr[N-1] && arr[0] <= arr[N/2])) pivotindex = 0;
	else if ((arr[N/2] >= arr[0] && arr[N/2] <= arr[N-1])
			|| (arr[N/2] >= arr[N-1] && arr[N/2] <= arr[0])) pivotindex = N/2;
	else pivotindex = N-1;
	pivotvalue=arr[pivotindex];

	// move pivot to end of array
	std::swap(arr[pivotindex],arr[N-1]);
	SwapPointsInArray(&pointarray[pivotindex],&pointarray[N-1]);
	newpivotindex=0;

	// partition list and array
	for(i=0;i<N;++i){
		if(arr[i] <= pivotvalue){
			std::swap(arr[newpivotindex],arr[i]);
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
	
	if(N <= 1) return ;

	// pick pivot as the median of the first, last and middle values
	if ((arr[0] >= arr[N/2] && arr[0] <= arr[N-1])
			|| (arr[0] >= arr[N-1] && arr[0] <= arr[N/2])) pivotindex = 0;
	else if ((arr[N/2] >= arr[0] && arr[N/2] <= arr[N-1])
			|| (arr[N/2] >= arr[N-1] && arr[N/2] <= arr[0])) pivotindex = N/2;
	else pivotindex = N-1;
	pivotvalue=arr[pivotindex];

	// move pivet to end of array
	std::swap(arr[pivotindex],arr[N-1]);
	std::swap(particles[pivotindex],particles[N-1]);
	newpivotindex=0;

	// partition list and array
	for(i=0;i<N;++i){
		if(arr[i] <= pivotvalue){
			std::swap(arr[newpivotindex],arr[i]);
			std::swap(particles[newpivotindex],particles[i]);
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

	*pivotindex=0;

	for(i=0;i<N;++i){
		if(arr[i] <= pivotvalue){
			std::swap(arr[*pivotindex],arr[i]);
			std::swap(particles[*pivotindex],particles[i]);
			++(*pivotindex);
		}
	}

	return ;
}
void quickPartitionPoints(double pivotvalue,unsigned long *pivotindex
		,Point *pointarray,double *arr,unsigned long N){
	unsigned long i;

	*pivotindex=0;

	for(i=0;i<N;++i){
		if(arr[i] <= pivotvalue){
			std::swap(arr[*pivotindex],arr[i]);
			SwapPointsInArray(&pointarray[*pivotindex],&pointarray[i]);
			++(*pivotindex);
		}
	}

	return ;
}

// return 1 (0) if box is (not) within rmax of ray
int cutbox(PosType *ray,PosType *p1,PosType *p2,float rmax){
	/*  returns:  0 if whole box is outside rmax from ray[]
	 *            1 if whole box is inside circle but ray is not in the box
	 *            2 if ray[] is inside box
	 *            3 if box intersects circle but ray[] is not inside box
	 */
  short i,tick=0;
  double close[2],rtmp;
  double tmp1,tmp2;

  // find closest point on box borders to ray[]
  for(i=0;i<2;++i){
    if( ray[i] < p1[i] ){
      close[i]=p1[i];
    }else if(ray[i] > p2[i]){
      close[i]=p2[i];
    }else{
      close[i]=ray[i];
      ++tick;
    }
  }

  if(tick==2) return 2;  // ray is inside box

  for(i=0,rtmp=0;i<2;++i) rtmp += pow(ray[i] - close[i],2);

  if(rtmp>rmax*rmax) return 0;  // box is all outside circle

  // find farthest point on box border from ray[]
  for(i=0,rtmp=0;i<2;++i) rtmp += ((tmp1 = pow(ray[i]-p1[i],2)) > (tmp2=pow(ray[i]-p2[i],2))) ? tmp1 : tmp2;
  //for(i=0,rtmp=0;i<2;++i) rtmp += DMAX(pow(ray[i]-p1[i],2),pow(ray[i]-p2[i],2));

  if(rtmp<rmax*rmax) return 1;  // box is all inside circle

  return 3;  // box intersects circle
}
}
