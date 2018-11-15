
/* double_sort is adapted from NR sort2 to take an unsigned long */
/* second array brr */
/* note #undef's at end of file */
#include "slsimlib.h"

#include <nrutil.h>
#include <mutex>
#include <thread>
#include <future>


#define NRANSI
#define M 7
#define NSTACK 50
namespace Utilities{
  void double_sort(unsigned long n, PosType *arr, unsigned long *brr)
  {
    unsigned long i,ir=n,j,k,l=1,*istack,b;
    int jstack=0;
    PosType a;
    
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
  void double_sort_points(unsigned long n, PosType *arr, Point *brr){
    unsigned long i,ir=n,j,k,l=1,*istack;
    long jstack=0;
    PosType a;
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
  
  void quicksortPoints(Point *pointarray,PosType *arr,unsigned long N){
    PosType pivotvalue;
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
    if(newpivotindex != 0) --newpivotindex;
    
    quicksortPoints(pointarray,arr,newpivotindex);
    quicksortPoints(&pointarray[newpivotindex+1],&arr[newpivotindex+1],N-newpivotindex-1);
    
    return ;
  }
  
  void quicksortPoints(Point *pointarray,double (*func)(Point &),unsigned long N){
    PosType pivotvalue;
    unsigned long pivotindex,newpivotindex,i;
    
    if(N <= 1) return ;
    
    // pick pivot as the median of the first, last and middle values
    if ((func(pointarray[0]) >= func(pointarray[N/2]) && func(pointarray[0]) <= func(pointarray[N-1]))
        || (func(pointarray[0]) >= func(pointarray[N-1]) && func(pointarray[0]) <= func(pointarray[N/2]))) pivotindex = 0;
    else if ((func(pointarray[N/2]) >= func(pointarray[0]) && func(pointarray[N/2]) <= func(pointarray[N-1]))
             || (func(pointarray[N/2]) >= func(pointarray[N-1]) && func(pointarray[N/2]) <= func(pointarray[0]))) pivotindex = N/2;
    else pivotindex = N-1;
    pivotvalue=func(pointarray[pivotindex]);
    
    // move pivot to end of array
    SwapPointsInArray(&pointarray[pivotindex],&pointarray[N-1]);
    newpivotindex=0;
    
    // partition list and array
    for(i=0;i<N;++i){
      if(func(pointarray[i]) <= pivotvalue){
        SwapPointsInArray(&pointarray[newpivotindex],&pointarray[i]);
        ++newpivotindex;
      }
    }
    if(newpivotindex != 0) --newpivotindex;
    
    quicksortPoints(pointarray,func,newpivotindex);
    quicksortPoints(&pointarray[newpivotindex+1],func,N-newpivotindex-1);
    
    return ;
  }

  void quicksort(unsigned long *particles,PosType *arr,unsigned long N){
    
    std::vector<size_t> index(N);
    
    Utilities::sort_indexes(arr,index,N);
    
    Utilities::apply_permutation(particles,index);
    Utilities::apply_permutation(arr,index);
    
    //std::cout << arr[0] << " " << arr[1] << " " << arr[2] << std::endl;
    assert(arr[0] <= arr[N-1]);
    return;

    
    PosType pivotvalue;
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
    if(newpivotindex != 0) --newpivotindex;
    
    quicksort(particles,arr,newpivotindex);
    quicksort(&particles[newpivotindex+1],&arr[newpivotindex+1],N-newpivotindex-1);
    
    return ;
  }
  
  /*
   * Partitions arr[] and particles[] into those with x <= pivotvalue and those with
   *  x > pivotvalue  pivotindex is left at first array value with x > pivotvalue
   */
  void quickPartition(PosType pivotvalue,unsigned long *pivotindex,unsigned long *particles
                      ,PosType *arr,unsigned long N){
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
  void quickPartitionPoints(PosType pivotvalue,unsigned long *pivotindex
                            ,Point *pointarray,PosType *arr,unsigned long N){
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
  void quickPartitionPoints(PosType pivotvalue,unsigned long *pivotindex
                            ,Point *pointarray,PosType (*func)(Point &p),unsigned long N){
    unsigned long i;
    
    *pivotindex=0;
    
    for(i=0;i<N;++i){
      if(func(pointarray[i]) <= pivotvalue){
        SwapPointsInArray(&pointarray[*pivotindex],&pointarray[i]);
        ++(*pivotindex);
      }
    }
    
    return ;
  }
 
}
