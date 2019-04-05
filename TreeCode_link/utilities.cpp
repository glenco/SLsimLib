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

  s_points = NewPointArray(Npoints);

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

void log_polar_grid(Point *i_points,PosType rmax,PosType rmin,PosType *center,long Ngrid){
  long i;
  PosType r,theta;
  static long id=0;

  for(i=0;i<Ngrid*Ngrid;++i){
	  theta= 2*PI*(i % Ngrid)/Ngrid;
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
long IndexFromPosition(PosType *x,long Npixels,PosType range,const PosType *center){
	long ix,iy;
	PosType fx, fy;

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
  
  /// this is the nonsquare version of the function, it will return -1 is outside of region
  long IndexFromPosition(PosType *x,long Nx,long Ny,PosType range,const PosType *center){
    long ix,iy;
    PosType fx, fy;
    
    fx = (x[0] - center[0])/range;
    fy = (x[1] - center[1])*(Nx-1)/range/(Ny-1);
    if(fabs(fx) > 0.5) return -1;
    if(fabs(fy) > 0.5) return -1;
    
	  fx = (fx + 0.5)*(Nx-1) + 0.5;
	  fy = (fy + 0.5)*(Ny-1) + 0.5;
    
	  if (fx < 0.) ix = -1;
	  else ix = (long)(fx);
    
	  if (fy < 0.) iy = -1;
	  else iy = (long)(fy);
    
	  if( (ix>-1)*(ix<Nx) && (iy>-1)*(iy<Ny) ) return ix+Nx*iy;
	  return -1;
  }

  /// This should work for square regions
  void PositionFromIndex(unsigned long i,PosType *x,long Npixels,PosType range,PosType const *center){
    if(Npixels == 1){
      x[0] = center[0];
      x[1] = center[1];
      return;
    }
    x[0] = center[0] + range*( 1.0*(i%Npixels)/(Npixels-1) - 0.5 );
    x[1] = center[1] + range*( 1.0*(i/Npixels)/(Npixels-1) - 0.5 );
    return;
  }
  /// This should work for square or rectangular regions as long as Npixels and range are the x-axis values and the pixels are square
  void PositionFromIndex(unsigned long i,PosType *x,long Nx,long Ny,PosType range,PosType const *center){
    if(Nx == 1){
      x[0] = center[0];
      x[1] = center[1];
      return;
    }
    x[0] = center[0] + range*( 1.0*(i%Nx)/(Nx-1) - 0.5 );
    x[1] = center[1] + range*( 1.0*(i/Nx)/(Nx-1) ) - 0.5*range*(Ny-1)/(Nx-1);
    return;
  }

// 1d version
long IndexFromPosition(PosType x,long Npixels,PosType range,PosType center){
	PosType fx;
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
   *  After it is used once, later calls can use TwoDInterpolator(PosType *map) for the same point 
   *  in the same coordinate system to save time in calculating the indexes.
   */
  PosType TwoDInterpolator(
                          PosType *x
                          ,int Npixels
                          ,PosType range
                          ,PosType *center
                          ,PosType *map
                          ,bool init
                          ){
    static PosType fx, fy;
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
  PosType TwoDInterpolator(PosType *map){
    PosType *dummy = NULL;
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


  PosType **PosTypeMatrix(size_t rows, size_t cols)
  {
    PosType **matrix = new PosType*[rows];
    matrix[0] = new PosType[rows*cols];
    for (long i = 1; i < rows; ++i)
      matrix[i] = matrix[0] + i * cols;
    
    return matrix;
  }
  
  
  void free_PosTypeMatrix(PosType **matrix, size_t rows, size_t cols){
    if(rows) delete[] matrix[0];
    delete[] matrix;
  }
  
  PosType **PosTypeMatrix(long rows1,long rows2, long cols1,long cols2)
  {
    PosType **matrix = new PosType*[(rows2-rows1)+1] - rows1;
    matrix[rows1] = new PosType[(rows2-rows1+1)*(cols2-cols1+1)];
    for (long i = 1; i < (rows2-rows1+1); ++i)
      matrix[rows1 + i] = matrix[rows1] + i * (cols2-cols1+1);
    
    return matrix;
  }
  
  
  void free_PosTypeMatrix(PosType **matrix, long rows1,long rows2, long cols1,long cols2){
    if(rows2-rows1) delete[] matrix[rows1];
    delete[] &matrix[rows1];
  }
  
  
  /// convert (x,y) to d
  int HilbertCurve::xy2d (int x, int y) {
    int rx, ry, s, d=0;
    for (s=n/2; s>0; s/=2) {
      rx = (x & s) > 0;
      ry = (y & s) > 0;
      d += s * s * ((3 * rx) ^ ry);
      rot(s, &x, &y, rx, ry);
    }
    return d;
  }
  
  ///convert d to (x,y)
  void HilbertCurve::d2xy(int d, int *x, int *y) {
    int rx, ry, s, t=d;
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
  int HilbertCurve::xy2d (PosType x, PosType y) {
    int rx, ry;
    
    rx = (int)((x-xo[0])*n/range + 0.5);
    ry = (int)((y-xo[1])*n/range + 0.5);
    
    if(rx > n || ry > n) throw std::runtime_error("Point out of bounds.");
    return xy2d(rx,ry);
  }
  
  ///convert d to (x,y)
  void HilbertCurve::d2xy(int d, PosType *x, PosType *y) {
    int rx, ry;
    
    if(d > n*n) throw std::runtime_error("Point out of bounds.");
      
    d2xy(d,&rx,&ry);
    
    *x = xo[0] + range*rx;
    *y = xo[1] + range*ry;
  }
  
  //rotate/flip a quadrant appropriately
  void HilbertCurve::rot(int s,int *x, int *y, int rx, int ry) {
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
  
#ifdef ENABLE_CLANG
  RandomNumbers::RandomNumbers(unsigned int seed)
  {        
    rand_gen = std::mt19937(seed);
  }
  
  RandomNumbers::~RandomNumbers(void){
    //delete rand_gen;
  }
  
  /// return a uniform random number between 0 and 1
  PosType RandomNumbers::operator()(void){
    return std::generate_canonical<PosType, 10>(rand_gen);
  }
#endif
  
  RandomNumbers_NR::RandomNumbers_NR(long seed):
  calls(0),IM1(2147483399),IM2(2147483399),IA1(40014),IA2(40692),IQ1(53668),
  IQ2(52774),IR1(12211),IR2(3791),EPS(1.2e-7),idum2(123456789),iy(0),
    count(true),firstseed(seed)
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
  PosType RandomNumbers_NR::operator()(void){
    ++calls;
    return ran2();
  }
  
  PosType RandomNumbers_NR::ran2(void){
    long j;
    long k;
    PosType temp;
    
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
  
  
#ifdef ENABLE_FFTW
  
  
  void powerspectrum2d(
                       std::valarray<double> const &aa
                       ,std::valarray<double> const &bb
                       ,int nx                       
                       ,int ny
                       ,double boxlx
                       ,double boxly
                       ,std::vector<double> &ll
                       ,std::vector<double> &Pl
                       ,double zeropaddingfactor
                       )
  {
    // go in the fourir space doing the zero padding
    int zerosize = 4;
    
    int nl = ll.size();
    Pl.resize(nl);
    
    // size of the new map in x and y directions, factor by which each size is increased
    int Nnx=int(zerosize*nx);
    int Nny=int(zerosize*ny);
    double Nboxlx = boxlx*zerosize;
    double Nboxly = boxly*zerosize;
    
    std:: valarray<float> Na,Nb;
    Na.resize( Nnx*Nny );
    Nb.resize( Nnx*Nny );
    
    // assume locate in a rectangular map and build up the new one
    for( int i=0; i<Nnx; i++ ) for( int j=0; j<Nny; j++ ){
      Na[i+Nnx*j]=0;
      Nb[i+Nnx*j]=0;
      if(i>=int(Nnx/2-nx/2) && i<int(Nnx/2+nx/2) && j>=int(Nny/2-ny/2) && j<int(Nny/2+ny/2)){
        int ii = i-int(Nnx/2-nx/2);
        int jj = j-int(Nny/2-ny/2);
        
        if(ii>=nx || jj>=ny){
          std::cout << " 1 error mapping " << ii << "  " << jj << std::endl;
          exit(1);
        }
        if(ii<0 || jj<0){
          std::cout << " 2 error mapping " << ii << "  " << jj << std::endl;
          exit(1);
        }
        Na[i+Nnx*j]=aa[ii+nx*jj];
        Nb[i+Nnx*j]=bb[ii+nx*jj];
      }
    }
    
    // estimate the fourier transform
    // now we go in the Fourier Space
    double *dNa=new double[Nnx*Nny];
    double *dNb=new double[Nnx*Nny];
    
    fftw_complex *fNa=new fftw_complex[Nny*(Nnx/2+1)];
    fftw_complex *fNb=new fftw_complex[Nny*(Nnx/2+1)];
    fftw_complex *Nfcc=new fftw_complex[Nny*(Nnx/2+1)];
    //size_t *ik=new size_t[];
    fftw_complex *output=new fftw_complex[Nny*(Nnx/2+1)];
    for(int i=0;i<Nnx;i++) for(int j=0;j<Nny; j++){
      dNa[i+Nnx*j] = double(Na[i+Nnx*j]);
      dNb[i+Nnx*j] = double(Nb[i+Nnx*j]);
    }
    
    fftw_plan p1 = fftw_plan_dft_r2c_2d(Nny,Nnx,dNa,fNa,FFTW_ESTIMATE);
    fftw_execute( p1 );
    fftw_destroy_plan(p1);
    fftw_plan p2 = fftw_plan_dft_r2c_2d(Nny,Nnx,dNb,fNb,FFTW_ESTIMATE);
    fftw_execute( p2 );
    fftw_destroy_plan(p2);
    
    //double *ks=new double[Nny*(Nnx/2+1)];
    std::vector<double> ks(Nny*(Nnx/2+1));
    
    // fb and fa are then used to compute the power spectrum
    // build modes
    for( int i=0; i<Nnx/2+1; i++ ){
      // kx = i if i<n/2 else i-n
      // double kx=(i<Nnx/2)?double(i):double(i-Nnx);
      double kx=double(i);
      for( int j=0; j<Nny; j++ ){
        // double ky=double(j);
        double ky=(j<Nny/2)?double(j):double(j-Nny);
        // rescale respect to the box size
        ks[i+(Nnx/2+1)*j] = sqrt(kx*kx/Nboxlx/Nboxlx + ky*ky/Nboxly/Nboxly)*2.*M_PI;
        Nfcc[i+(Nnx/2+1)*j][0] =  (fNa[i+(Nnx/2+1)*j][0]*fNb[i+(Nnx/2+1)*j][0] +
                                   fNa[i+(Nnx/2+1)*j][1]*fNb[i+(Nnx/2+1)*j][1])*pow(1./2/M_PI/Nnx/Nny,2)*Nboxlx*Nboxly;
        Nfcc[i+(Nnx/2+1)*j][1] = -(fNa[i+(Nnx/2+1)*j][0]*fNb[i+(Nnx/2+1)*j][1] -
                                   fNa[i+(Nnx/2+1)*j][1]*fNb[i+(Nnx/2+1)*j][0])*pow(1./2/M_PI/Nnx/Nny,2)*Nboxlx*Nboxly;
      }
    }
    std::vector<size_t> ik(Nny*(Nnx/2+1));
    //gsl_sort_index (ik,ks,1,Nny*(Nnx/2+1));
    Utilities::sort_indexes<double>(ks,ik);
    
    double *nk=new double[Nny*(Nnx/2+1)];
    double *nc=new double[Nny*(Nnx/2+1)];
    // sorted vectors
    for(int i=0; i<Nnx/2+1; i++)
      for(int j=0; j<Nny;j++){
        nk[i+(Nnx/2+1)*j] = ks[ik[i+(Nnx/2+1)*j]];
        nc[i+(Nnx/2+1)*j] = Nfcc[ik[i+(Nnx/2+1)*j]][0];
      }
    std:: vector<double> bink(nl);
    // build the binned power spectrum
    Utilities::fill_linear(bink,nl,log10(nk[1]),log10(nk[Nny*(Nnx/2+1)-1]));
    double lk1,lk2;
    for(int i=0;i<nl;i++){
      if(i==0) lk1=bink[i];
      else lk1=bink[i]-0.5*(bink[i]-bink[i-1]);
      if(i==nl-1) lk2=bink[i];
      else lk2=bink[i]+0.5*(bink[i+1]-bink[i]);
      Pl[i]=0.;
      ll[i]=0.;
      int nin=0;
      // start from 1 because the first is 0
      for(int j=1;j<Nny*(Nnx/2+1)-1;j++){
        if(log10(nk[j])>=lk1 && log10(nk[j])<lk2){
          Pl[i]=Pl[i]+nc[j];
          ll[i]=ll[i]+nk[j];
          nin=nin+1;
        }
      }
      if(nin>0){
        Pl[i]=Pl[i]/double(nin)*zerosize*zerosize;
        ll[i]=ll[i]/double(nin);
      }
    }
    delete [] dNa;
    delete [] dNb;
    delete [] output;
    delete [] fNa;
    delete [] fNb;
    delete [] Nfcc;
    delete [] nk;
    delete [] nc;
  }
  
#endif
  
  std::valarray<double> AdaptiveSmooth(const std::valarray<double> &map_in,size_t Nx,size_t Ny,double value){
    
    std::valarray<double> map_out(map_in.size());
    long r,area;
    double val;
    for(long i=0;i<Nx;++i){
      for(long j=0;j<Ny;++j){
        r = 0;
        val = map_in[i+j*Nx];
        while(val < value && r < std::min(Nx, Ny) ){
          
          area = 0;
          val = 0;
          long imin,imax,jmin,jmax;
          
          imin = (i-r < 0) ? 0 : i-r;
          imax = (i+r > Nx-1) ? Nx-1 : i+r;
          
          jmin = (j-r < 0) ? 0 : j-r;
          jmax = (j+r > Ny-1) ? Ny-1 : j+r;
          
          
          for(long ii=imin;ii<=imax;++ii){
            for(long jj=jmin;jj<=jmax;++jj){
              if( (ii-i)*(ii-i) + (jj-j)*(jj-j) < r){
                val += map_in[ii+jj*Nx];
                ++area;
              }
            }
          }
          ++r;
        }
        
        map_out[i+j*Nx] = val/area;
      }
    }
    
    return map_out;
  }
  /** \brief Smooth a 2 dimensional map stored in a valarray with a density dependent kernel.
   
   The smoothing is done by finding the circle around each point whose total pixel values are larger than value.  In the case of a density map made from particles if value = (mass of particle)*(number of neighbours) an approximate N nearest neighbour smoothing is done.
   The
   **/
  std::vector<double> AdaptiveSmooth(const std::vector<double> &map_in,size_t Nx,size_t Ny,double value){
    
    std::vector<double> map_out(map_in.size());
    long r,area;
    double val;
    for(long i=0;i<Nx;++i){
      for(long j=0;j<Ny;++j){
        r = 0;
        val = map_in[i+j*Nx];
        while(val < value && r < std::min(Nx, Ny) ){
          
          area = 0;
          val = 0;
          long imin,imax,jmin,jmax;
          
          imin = (i-r < 0) ? 0 : i-r;
          imax = (i+r > Nx-1) ? Nx-1 : i+r;
          
          jmin = (j-r < 0) ? 0 : j-r;
          jmax = (j+r > Ny-1) ? Ny-1 : j+r;
          
          
          for(long ii=imin;ii<=imax;++ii){
            for(long jj=jmin;jj<=jmax;++jj){
              if( (ii-i)*(ii-i) + (jj-j)*(jj-j) < r){
                val += map_in[ii+jj*Nx];
                ++area;
              }
            }
          }
          ++r;
        }
        
        map_out[i+j*Nx] = val/area;
      }
    }
    
    return map_out;
  }

  int GetNThreads(){return N_THREADS;}
}

void Utilities::IO::ReadFileNames(
                   std::string dir              /// path to directory containing fits files
                   ,const std::string filespec /// string of charactors in file name that are matched. It can be an empty string.
                   ,std::vector<std::string> & filenames  /// output vector of PixelMaps
                   ,bool verbose){
  
  DIR *dp = opendir( dir.c_str() );
  struct dirent *dirp;
  struct stat filestat;
  std::string filepath,filename;
  
  if (dp == NULL)
  {
    std::cerr << "Cannot find directory" << std::endl;
    throw std::runtime_error("error opening directory");
    return;
  }
  
  while ((dirp = readdir( dp )) )
  {
    filepath = dir + "/" + dirp->d_name;
    
    // If the file is a directory (or is in some way invalid) we'll skip it
    if (stat( filepath.c_str(), &filestat )) continue;
    if (S_ISDIR( filestat.st_mode ))         continue;
    
    filename = dirp->d_name;
    if(filename.find(filespec) !=  std::string::npos){
      if(verbose) std::cout << "adding " << filepath << std::endl;
      filenames.push_back(filename);
    }
  }
  
  closedir( dp );
  
  std::cout << filenames.size() << " file names." << std::endl;
  return ;
}

/// Count the number of columns in a ASCII data file
int Utilities::IO::CountColumns(std::string filename,char comment_char
                                ,char deliniator
){
  
  std::ifstream file(filename);
  // find number of particles
  if (!file.is_open()){
    std::cerr << "file " << filename << " cann't be opened." << std::endl;
    throw std::runtime_error("no file");
  }
  
  std::string line;
  // read comment lines and first data line
  do{
    std::getline(file,line);
    if(!file) break;  // probably EOF
  }while(line[0] == comment_char);
  
  return  NumberOfEntries(line,deliniator);
}

int Utilities::IO::NumberOfEntries(const std::string &string,char deliniator){
  size_t number = 0;
  auto it = string.begin();
  while(it != string.end()){
    while( it != string.end() && (*it == ' ' || *it == deliniator)) ++it;
    if(it != string.end()) ++number;
    while(it != string.end() && *it != deliniator) ++it;
  }
  return number;
}

Utilities::XYcsvLookUp::XYcsvLookUp(
                                    std::string datafile   /// input catalog file in csv format
                                    ,std::string Xlabel
                                    ,std::string Ylabel
                                    ,int Nxbins  /// number of X bins
                                    ,size_t MaxNumber
                                    ,bool verbose)
:filename(datafile)
{
  Utilities::IO::ReadCSVnumerical2<double>(datafile,data,column_names,MaxNumber);
  if(verbose){
    for(auto name : column_names){ std::cout << name << " " ;}
    std::cout << std::endl;
    for(int i=0 ; i < 1 ; ++i){
      for(auto v : data[i] ) std::cout << v << " " ;
      std::cout << std::endl;
    }
    std::cout << data.size() << " rows with " << column_names.size() << " columns read." << std::endl;
  }
  
  int i=0;
  Xindex = Yindex = -1;
  for(auto name : column_names){
    if(name == Xlabel) Xindex = i;
    if(name == Ylabel) Yindex = i;
    ++i;
  }
  if(Xindex == -1){
    std::cerr << filename << " needs a column named " << Xlabel << std::endl;
    throw std::invalid_argument("No column named: " + Xlabel);
  }
  if(Yindex == -1){
    std::cerr << filename << " needs a column named " << Ylabel << std::endl
    << " They are :" << std::endl;
    for(auto c : column_names ) std::cout << c << " ";
    std::cout << std::endl;
    throw std::invalid_argument("No column named: " + Ylabel);
  }
  
  //************* test ********************
  //xmin = 100;
  //xmax = -1;
  //for(auto d : data){
  //  if(d[Xindex] < xmin) xmin = d[Xindex];
  //  if(d[Xindex] > xmax) xmax = d[Xindex];
  //}
  //**********************************
  
  // sort by redshift
  std::sort(data.begin(),data.end(), [this](const std::vector<double> &v1,const std::vector<double> &v2){return v1[Xindex] < v2[Xindex];});
  //&v1,std::vector<double> &v2){return v1[1] < v2[1];});
  NinXbins = data.size()/Nxbins;
  Xborders.resize(Nxbins);
  Xborders[0]=0.0;
  xmin =  data[0][Xindex];
  xmax =  data.back()[Xindex];
  if(verbose){
    std::cout << "min X : "<< data[0][Xindex] << " max X : "
    << data.back()[Xindex] << std::endl;
    std::cout << column_names[Xindex] << " Bins " << std::endl;
  }
  for(int i=1 ; i<Nxbins ; ++i){
    Xborders[i] = data[i*NinXbins][Xindex];
    if(verbose) std::cout << i << " " << Xborders[i] << std::endl;
  }
  
  // set up iterators to boundaries of x bins
  borders.resize(Nxbins + 1);
  borders[0] = data.begin();
  borders.back() = data.end();
  for(int i=1 ; i < Nxbins ; ++i) borders[i] = borders[i-1] + NinXbins;
  
  // sort by mass within x bins
  for(int i=0 ; i < Nxbins ; ++i){
    std::sort(borders[i],borders[i+1], [this](const std::vector<double> &v1
                                              ,const std::vector<double> &v2){return v1[Yindex] < v2[Yindex];});
  }
  
  current = data.begin();
}

Utilities::XYcsvLookUp::XYcsvLookUp(
                                    std::string datafile   /// input catalog file in csv format
                                    ,std::string Xlabel
                                    ,std::string Ylabel
                                    ,std::vector<double> Xbins
                                    ,size_t MaxNumber
                                    ,bool verbose)
:Xborders(Xbins),filename(datafile)
{
  
  Utilities::IO::ReadCSVnumerical2<double>(datafile,data,column_names,MaxNumber);
  if(verbose){
    for(auto name : column_names){ std::cout << name << " " ;}
    std::cout << std::endl;
    for(int i=0 ; i < 1 ; ++i){
      for(auto v : data[i] ) std::cout << v << " " ;
      std::cout << std::endl;
    }
    std::cout << data.size() << " rows with " << column_names.size() << " columns read." << std::endl;
  }
  
  int i=0;
  Xindex = Yindex = -1;
  for(auto name : column_names){
    if(name == Xlabel) Xindex = i;
    if(name == Ylabel) Yindex = i;
    ++i;
  }
  if(Xindex == -1){
    std::cerr << filename << " needs a column named " << Xlabel << std::endl;
    throw std::invalid_argument("No column named: " + Xlabel);
  }
  if(Yindex == -1){
    std::cerr << filename << " needs a column named " << Ylabel << std::endl
    << " They are :" << std::endl;
    for(auto c : column_names ) std::cout << c << " ";
    std::cout << std::endl;
    throw std::invalid_argument("No column named: " + Ylabel);
  }
  
  // sort by redshift
  std::sort(data.begin(),data.end(), [this](const std::vector<double> &v1
                                            ,const std::vector<double> &v2){return v1[Xindex] < v2[Xindex];});

  size_t Nxbins = Xborders.size();
  NinXbins = data.size()/Nxbins;
  xmin =  data[0][Xindex];
  xmax =  data.back()[Xindex];
  if(verbose){
    std::cout << "min X : "<< data[0][Xindex] << " max X"
    << data.back()[Xindex] << std::endl;
    std::cout << "redshift bins " << std::endl;
  }
  
  // set up iterators to boundaries of x bins
  borders.resize(Nxbins + 1);
  borders[0] = data.begin();
  borders.back() = data.end();
  for(int i=1 ; i < Nxbins ; ++i) borders[i] = borders[i-1] + NinXbins;
  
  // sort by mass within x bins
  for(int i=0 ; i < Nxbins ; ++i){
    std::sort(borders[i],borders[i+1], [this](const std::vector<double> &v1
                                              ,const std::vector<double> &v2){return v1[Yindex] < v2[Yindex];});
  }
  
  current = data.begin();
}

std::vector<double> Utilities::XYcsvLookUp::find(double x,double y){
  long xbin = Utilities::locate(Xborders, x);
  current = std::upper_bound(borders[xbin],borders[xbin+1],y
                             , [this](double y,const std::vector<double> &v1){return y < v1[Yindex];});
  
  if(current == borders[xbin+1]) current = borders[xbin+1] - 1;
  //std::cout << "XYcsvLookUp boundaries :" << (*borders[xbin])[Yindex] << " " << (*borders[xbin+1])[Yindex] << std::endl;
  //assert(fabs(y-(*current)[Yindex]) < 1);
  return *current;
}
double Utilities::XYcsvLookUp::Ymin(double x) const{
  long xbin = Utilities::locate(Xborders, x);
  return (*borders[xbin])[Yindex];
}
double Utilities::XYcsvLookUp::Ymax(double x) const{
  long xbin = Utilities::locate(Xborders, x);
  return (*(borders[xbin+1]-1))[Yindex];
}

double Utilities::XYcsvLookUp::operator[](std::string label) const{
  int i = 0;
  for(auto name : column_names){
    if(name == label){
      return (*current)[i];
    }
    ++i;
  }
  for(auto c : column_names ) std::cout << c << " ";
  std::cout << std::endl;
  throw std::invalid_argument(label + " was not one of the columns of the galaxy data file :" + filename);
}

void Utilities::splitstring(std::string &line,std::vector<std::string> &vec
                            ,const std::string &delimiter){
  size_t pos = 0;
  
  while ((pos = line.find(delimiter)) != std::string::npos) {
    vec.push_back(line.substr(0, pos));
    line.erase(0, pos + delimiter.length());
  }
  vec.push_back(line.substr(0, pos));
}
