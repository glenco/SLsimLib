/**
 * utilities.c
 *
 *  Created on: Sep 8, 2009
 *      Author: R.B. Metcalf
 */


#include <slsimlib.h>

#define DLTABLE	256
std::vector<double> zz, da;

Point *LinkToSourcePoints(Point *i_points,unsigned long Npoints){
  Point *s_points;
  long i;

  if(Npoints < 1) return NULL;

  s_points = NewPointArray(Npoints,true);

#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
  for(i=0;i<Npoints;++i){
    s_points[i].id=i_points[i].id;
      // link images and source points
    i_points[i].image=&s_points[i];
    s_points[i].image=&i_points[i];
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
 */

long IndexFromPosition(double *x,long Npixels,double range,double *center){
	long ix,iy;

	  ix=(long)( ((x[0] - center[0])/range + 0.5)*(Npixels-1)+0.5);
	  iy=(long)( ((x[1] - center[1])/range + 0.5)*(Npixels-1)+0.5);

/*	  std::printf("point %e %e  map %e %e\n",x[0],x[1]
      ,(ix*1.0/(Npixels-1.) - 0.5)*range + center[0]
      ,(iy*1.0/(Npixels-1.) - 0.5)*range + center[1]);
*/

	  if( (ix>-1)*(ix<Npixels) && (iy>-1)*(iy<Npixels) ) return ix+Npixels*iy;
	  return -1;
}
/** \ingroup Utill
 *
 */
void PositionFromIndex(unsigned long i,double *x,long Npixels,double range,double *center){
    x[0] = center[0] + range*( 1.0*(i%Npixels)/(Npixels-1) - 0.5 );
    x[1] = center[1] + range*( 1.0*(i/Npixels)/(Npixels-1) - 0.5 );
    return;
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


void makeDlTable(CosmoHndl cosmo, double zmax){
	fill_linear(zz,DLTABLE,0.1,zmax);

	da.resize(DLTABLE);

	for(int i = 0; i < DLTABLE; i++){
		da[i] = cosmo->angDist(0, zz[i]);
	}
}
