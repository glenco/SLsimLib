/**
 * utilities.c
 *
 *  Created on: Sep 8, 2009
 *      Author: R.B. Metcalf
 */


#include "slsimlib.h"
#include "utilities_slsim.h"
#include <random>

Point *LinkToSourcePoints(Point *i_points,unsigned long Npoints){
  Point *s_points;
  long i;

  if(Npoints < 1) return NULL;

  s_points = NewPointArray(Npoints,true);

  for(i=0;i<Npoints;++i){
    s_points[i].id=i_points[i].id;
      // link images and source points
    i_points[i].image=&s_points[i];
    s_points[i].image=&i_points[i];
    
    s_points[i].gridsize = i_points[i].gridsize;
  }

  return s_points;
}

// Calculates the area in the image by simply adding up all the cells in the imageinfo->point array
void findarea(OldImageInfo *imageinfo){
  unsigned long j;

  for(j=0,imageinfo->area=0,imageinfo->area_error=0;j<imageinfo->Npoints;++j){
	  if( fabs(imageinfo->points[j].leaf->boundary_p2[0] - imageinfo->points[j].leaf->boundary_p1[0]
	             - imageinfo->points[j].gridsize)/imageinfo->points[j].gridsize > 1.0e-4 ){
		  std::printf(" miss aligned gridsize %.5f %.5f \n",(imageinfo->points[j].leaf->boundary_p2[0] - imageinfo->points[j].leaf->boundary_p1[0])
		  	          /imageinfo->points[j].gridsize,(imageinfo->points[j].leaf->boundary_p2[1] - imageinfo->points[j].leaf->boundary_p1[1])
		  	          /imageinfo->points[j].gridsize);
	  }

    imageinfo->area += pow(imageinfo->points[j].gridsize,2);
    if(imageinfo->area_error < pow(imageinfo->points[j].gridsize,2))
      imageinfo->area_error = pow(imageinfo->points[j].gridsize,2);
  }

/*   MoveToTopList(imageinfo->outerborder); */
/*   for(j=0,imageinfo->area_error=0;j<imageinfo->outerborder->Npoints;++j){ */
/*     imageinfo->area_error += pow(imageinfo->outerborder->current->gridsize,2); */
/*     MoveDownList(imageinfo->outerborder); */
/*   } */
/*   imageinfo->area += 0.5*imageinfo->area_error; */

  imageinfo->area_error /= imageinfo->area;
}

/**  /ingroup Utill
 *
 *  make a new 2d grid with logarithmically distributed radii
 * Ngrid - both the radial and asmuthal number of points,
 * 	     total number of points is Ngrid*Ngrid
 * WARRNING: This should not be used as initial grid for adaptive
 *           grid calculations!  These need to be done on a rectolinear
 *           grid.  For testing the ray shooter only.
 */

void log_polar_grid(Point *i_points,double rmax,double rmin,double *center,long Ngrid){
  long i;
  double r,theta;
  static long id=0;

  for(i=0;i<Ngrid*Ngrid;++i){
	  theta= 2*pi*(i % Ngrid)/Ngrid;
	  r= rmin*exp( (i/Ngrid)*1.0/(Ngrid-1) * log(rmax/rmin) );
	  i_points[i].id=id;
      ++id;
      i_points[i].x[0] = center[0] + r*cos(theta);
      i_points[i].x[1] = center[1] + r*sin(theta);
      i_points[i].gridsize=1.0;
      //std::printf("x = %e %e\n",i_points[i].x[0],i_points[i].x[1]);
  }

  return;
}

/** \ingroup Utill
 *
 * The two functions below are inverses of each other for converting
 *   between a 1d array index and a square grid of positions
 *   Npixels in the number of point is 1 dimension
 *   index is between 0 and Npixels*Npixels-1
 *   If x is outside of the region -1 is returned.
 *
 */
namespace Utilities{
long IndexFromPosition(double *x,long Npixels,double range,double *center){
	long ix,iy;
	double fx, fy;

	  fx = ((x[0] - center[0])/range + 0.5)*(Npixels-1)+0.5;
	  fy = ((x[1] - center[1])/range + 0.5)*(Npixels-1)+0.5;

/*	  std::printf("point %e %e  map %e %e\n",x[0],x[1]
      ,(ix*1.0/(Npixels-1.) - 0.5)*range + center[0]
      ,(iy*1.0/(Npixels-1.) - 0.5)*range + center[1]);
*/
	  if (fx < 0.) ix = -1;
	  else ix = (long)(fx);

	  if (fy < 0.) iy = -1;
	  else iy = (long)(fy);

	  if( (ix>-1)*(ix<Npixels) && (iy>-1)*(iy<Npixels) ) return ix+Npixels*iy;
	  return -1;
}
  
/** \ingroup Utill
 *
 */
void PositionFromIndex(unsigned long i,double *x,long Npixels,double range,const double *center){
  if(Npixels == 1){
    x[0] = center[0];
    x[1] = center[1];
    return;
  }
  x[0] = center[0] + range*( 1.0*(i%Npixels)/(Npixels-1) - 0.5 );
  x[1] = center[1] + range*( 1.0*(i/Npixels)/(Npixels-1) - 0.5 );
  return;
}

// 1d version
long IndexFromPosition(double x,long Npixels,double range,double center){
	double fx;
	long ix;

	  fx = ((x - center)/range + 0.5)*(Npixels-1)+0.5;
	  if (fx < 0.) ix = -1;
	  else ix = (long)(fx);

	  if( (ix>-1)*(ix<Npixels)) return ix;
	  return -1;
}

  /** \ingroup Utill
   * \brief bilinear interpolation from a map.
   *
   *  Out of bounds points return 0.  map is a i dimensional array representing a 2 dimensional map.
   *  Don't use init.
   *  After it is used once, later calls can use TwoDInterpolator(double *map) for the same point 
   *  in the same coordinate system to save time in calculating the idndexes.
   */
  double TwoDInterpolator(
                          double *x
                          ,int Npixels
                          ,double range
                          ,double *center
                          ,double *map
                          ,bool init
                          ){
    static double fx, fy;
    static unsigned long index;
    static bool initialized = false;
    
    if(init){
      unsigned long ix,iy;
      fx = ((x[0] - center[0])/range + 0.5)*(Npixels-1);
      fy = ((x[1] - center[1])/range + 0.5)*(Npixels-1);
    
      if (fx < 0. || fx > Npixels-1) return 0.0;
      else ix = (unsigned long)(fx);
      if(ix == Npixels-1) ix = Npixels-2;
    
      if (fy < 0. || fx > Npixels-1) return 0.0;
      else iy = (unsigned long)(fy);
      if(iy == Npixels-1) iy = Npixels-2;
    
      index = ix + Npixels*iy;
      
      /** bilinear interpolation */
      fx = center[0] + range*( 1.0*(ix)/(Npixels-1) - 0.5 );
      fy = center[1] + range*( 1.0*(iy)/(Npixels-1) - 0.5 );
      fx=(x[0] - fx)*(Npixels-1);
      fy=(x[1] - fy)*(Npixels-1);

      initialized = true;
    }
    
    if(!initialized){
      ERROR_MESSAGE();
      std::cout << " TwoDInterpolator() needs to be initialized. " << std::endl;
      throw std::runtime_error("Not initialized");
    }
    
    return (1-fx)*(1-fy)*map[index] + fx*(1-fy)*map[index+Npixels] + fx*fy*map[index+1+Npixels]
           + (1-fx)*fy*map[index+1];
  }
  double TwoDInterpolator(double *map){
    double *dummy = NULL;
    return TwoDInterpolator(dummy,0,0,dummy,map,false);
  }


/** \ingroup Utill
 * This function finds the largest power of 2 that is < k
 */
unsigned long prevpower(unsigned long k){
	int i;
	if (k == 0)
		return 1;
	k--;
	for (i=1; i<sizeof(unsigned long)*sizeof(char); i<<=1)
		k = k | k >> i;
	return (k+1)/2;
}


PosType **PosTypeMatrix(long rows, long cols)
{
  PosType **matrix = new PosType*[rows];
  matrix[0] = new PosType[rows*cols];
  for (long i = 1; i < rows; ++i)
    matrix[i] = matrix[0] + i * cols;

  return matrix;
}


  void free_PosTypeMatrix(PosType **matrix, long rows, long cols){
    if(rows) delete[] matrix[0];
    delete[] matrix;
  }
  
  
  /// convert (x,y) to d
  long HilbertCurve::xy2d (long x, long y) {
    long rx, ry, s, d=0;
    for (s=n/2; s>0; s/=2) {
      rx = (x & s) > 0;
      ry = (y & s) > 0;
      d += s * s * ((3 * rx) ^ ry);
      rot(s, &x, &y, rx, ry);
    }
    return d;
  }
  
  ///convert d to (x,y)
  void HilbertCurve::d2xy(long d, long *x, long *y) {
    long rx, ry, s, t=d;
    *x = *y = 0;
    for (s=1; s<n; s*=2) {
      rx = 1 & (t/2);
      ry = 1 & (t ^ rx);
      rot(s, x, y, rx, ry);
      *x += s * rx;
      *y += s * ry;
      t /= 4;
    }
  }
  /// convert (x,y) to d
  long HilbertCurve::xy2d (double x, double y) {
    long rx, ry;
    
    rx = (long)((x-xo[0])*n/range + 0.5);
    ry = (long)((y-xo[1])*n/range + 0.5);
    
    if(rx > n || ry > n) throw std::runtime_error("Point out of bounds.");
    return xy2d(rx,ry);
  }
  
  ///convert d to (x,y)
  void HilbertCurve::d2xy(long d, double *x, double *y) {
    long rx, ry;
    
    if(d > n*n) throw std::runtime_error("Point out of bounds.");
      
    d2xy(d,&rx,&ry);
    
    *x = xo[0] + range*rx;
    *y = xo[1] + range*ry;
  }
  
  //rotate/flip a quadrant appropriately
  void HilbertCurve::rot(long s,long *x, long *y, long rx, long ry) {
    if (ry == 0) {
      if (rx == 1) {
        *x = s-1 - *x;
        *y = s-1 - *y;
      }
      
      //Swap x and y
      int t  = *x;
      *x = *y;
      *y = t;
    }
  }
  
  RandomNumbers::RandomNumbers(unsigned int seed)
  {        
    rand_gen = std::mt19937(seed);
  }
  
  RandomNumbers::~RandomNumbers(void){
    //delete rand_gen;
  }
  
  /// return a uniform random number between 0 and 1
  double RandomNumbers::operator()(void){
    return std::generate_canonical<double, 10>(rand_gen);
  }
  
  
  RandomNumbers_NR::RandomNumbers_NR(long seed):
  IM1(2147483399),IM2(2147483399),IA1(40014),IA2(40692),IQ1(53668),
  IQ2(52774),IR1(12211),IR2(3791),EPS(1.2e-7),idum2(123456789),iy(0)
  {
    
    AM = (1.0/IM1);
    //IMM1 = (IM1-1);
    NDIV = (1+(IM1-1)/32);
    RNMX = (1.0-EPS);
    
    if(seed > 0) seed *= -1;
    idum = seed;
    
    long k,j;
    
    if (-(idum) < 1) idum=1;
    else idum = -(idum);
    idum2 = idum;
    for (j=32+7;j>=0;j--) {
      k = idum/IQ1;
      idum=IA1*(idum-k*IQ1)-k*IR1;
      if (idum < 0) idum += IM1;
      if (j < 32) iv[j] = idum;
    }
    iy=iv[0];
    
  }
  
  /// return a uniform random number between 0 and 1
  double RandomNumbers_NR::operator()(void){
    return ran2();
  }
  
  double RandomNumbers_NR::ran2(void){
    long j;
    long k;
    double temp;
    
    k=(idum)/IQ1;
    idum=IA1*(idum-k*IQ1)-k*IR1;
    if (idum < 0) idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = idum;
    if (iy < 1) iy += IM1-1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
  }

}
