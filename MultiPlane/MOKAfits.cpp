/*
 * fits.cpp
 *
 *  Created on: Jun 19, 2012
 *      Author: mpetkova
 */

#include "MOKAlens.h"
#include <fstream>

#include "cpfits.h"

#ifdef ENABLE_FFTW
#include "fftw3.h"
#endif

// carlo test begin
/*
 int getiMIN(int a, int b){
 return (a) < (b) ? (a):(b);
 }
 
 int getiMAX(int a, int b){
 return (a) > (b) ? (a):(b);
 }
 
 double getMapVal(std::valarray<double> &map,int npix,double x,double y, std:: vector<double> &xi,std:: vector<double> &yi){
 size_t i=locate (xi,x);
 i=getiMIN( getiMAX( i, 0 ), npix-2 );
 size_t j=locate (yi,y);
 j=getiMIN( getiMAX( j, 0 ), npix-2 );
 size_t tnpix = npix;
 double fQ11 = map[i+tnpix*j];
 double fQ21;
 if(i<npix-1) fQ21 = map[(i+1)+tnpix*j];
 else fQ21 = 0;
 double fQ12;
 if(j<npix-1) fQ12 = map[i+tnpix*(j+1)];
 else fQ12 = 0;
 double fQ22;
 if(i<npix-1 && j<npix-1) fQ22 = map[(i+1)+tnpix*(j+1)];
 else fQ22 = 0;
 double fxy;
 //if(i<npix-1 && j<npix-1){
 fxy = fQ11*(xi[i+1]-x)*(yi[j+1]-y)/((xi[i+1]-xi[i])*(yi[j+1]-yi[j])) +
 fQ21*(x-xi[i])*(yi[j+1]-y)/((xi[i+1]-xi[i])*(yi[j+1]-yi[j])) +
 fQ12*(xi[i+1]-x)*(y-yi[j])/((xi[i+1]-xi[i])*(yi[j+1]-yi[j])) +
 fQ22*(x-xi[i])*(y-yi[j])/((xi[i+1]-xi[i])*(yi[j+1]-yi[j]));
 //}
 
 
 if(fxy!=fxy){
 std:: cout << x << "  " << y  << std:: endl;
 std:: cout << i << "  " << j << "  " << npix << std:: endl;
 std:: cout << "x = " << xi[i] << "  " << xi[i+1] << std:: endl;
 std:: cout << "y = " << yi[j] << "  " << yi[j+1] << std:: endl;
 std:: cout << fQ11 << "  " << fQ21 << "  " << fQ12 << "  " << fQ22 << std:: endl;
 }
 return fxy;
 }
 */
// carlo test end

std::vector<long> LensHaloMassMap::getDims(){
  std::vector<long> size(2);
  size[0] = map.nx;
  size[1] = map.ny;
  
  return size;
}

void LensHaloMassMap::writeImage(std::string filename){
  map.write(filename);
}

/**
 * routine used by fof to link nearby grid cell points
 */
void make_friendship(int ii,int ji,int np,std:: vector<int> &friends, std:: vector<double> &pointdist){
  for(int jj=0;jj<np;jj++){
    if(friends[ji+np*jj]!=0){
      if(friends[ji+np*jj]<0){
        friends[ii+np*jj]=-(ii+1);
      }
      else{
        friends[ii+np*jj]=(ii+1);
      }
      friends[ji+np*jj]=0;
    }
  }
  friends[ii+np*ji]=-(ii+1);
}

/*
 * given a a set of grid points xci and yci and an interpixeld distance l return the id of the
 * fof group nearest to the centre
 */

int fof(double l,std:: vector<double> xci, std:: vector<double> yci, std:: vector<int> &groupid){
  int np = xci.size();
  std:: vector<int> friends(np*np);
  std:: vector<double> pointdist(np*np);
  for(int ii = 0;ii<np; ii++) for(int ji = 0;ji<np; ji++){
    pointdist[ii+np*ji] = sqrt( pow(xci[ii] - xci[ji],2) + pow(yci[ii] - yci[ji],2));
    groupid[ii] = 0;
    friends[ii+np*ji]=0;
  }
  for(int ii=0;ii<np;ii++) for(int ji = 0;ji<np; ji++){
    if(pointdist[ii+np*ji]<=1.5*l) friends[ii+np*ji] = ii+1;
  }
  for(int ii=0;ii<np;ii++){
    int r = 0;
    while(r==0){
      r=1;
      for(int ji=0;ji<np;ji++){
        if(friends[ii+np*ji]>0){
          if(ii!=ji){
            make_friendship(ii,ji,np,friends,pointdist);
            r=0;
          }
        }
      }
    }
  }
  for(int ii=0;ii<np;ii++){
    int p=0;
    for(int ji=0;ji<np;ji++){
      if(friends[ji+np*ii]!=0) p++;
      if(p==2){
        std:: cout << ji << "  " << ii << ":  "  << friends[ji+np*ii] << "  " << friends[ii+np*ji] << std:: endl;
        exit(1);
      }
    }
  }
  // count the particles in each group
  int kt = 0;
  int ng= 0;
  std:: vector<double> distcentre;
  std:: vector<int> idgroup;
  for(int ii=0;ii<np;ii++){
    int k = 0;
    double xcm=0;
    double ycm=0;
    for(int ji=0;ji<np;ji++){
      if(friends[ii+np*ji]!=0){
        k++;
        groupid[ji]=ii+1;
        xcm += xci[ji];
        ycm += yci[ji];
      }
    }
    if(k>4){
      ng++;
      xcm/=k;
      ycm/=k;
      distcentre.push_back(sqrt(xcm*xcm+ycm*ycm));
      idgroup.push_back(ii+1);
      // std:: cout << "  " << ii+1 << "  " << k << "   " << sqrt(xcm*xcm+ycm*ycm) << std:: endl;
    }
    kt = kt + k;
  }
  if(kt != np){
    std:: cout << " number of screaned particles : " << kt << std:: endl;
    std:: cout << " differes from the number of particles : " << np << std:: endl;
    std:: cout << " number of group found : " << ng << std:: endl;
    std:: cout << "     " << std:: endl;
    std:: cout << " I will STOP here!!! " << std:: endl;
    exit(1);
  }
  if(idgroup.size()>0){
    std:: vector<double>::iterator it = min_element(distcentre.begin(), distcentre.end());
    int minpos = idgroup[distance(distcentre.begin(), it)];
    // std:: cout << " nearest to the centre " << minpos << std:: endl;
    return minpos;
  }
  else return 0;
  /* Make a histogram of the data */
  /*
   std::vector< int > histogram(np,0);
   std::vector< int >::iterator it = groupid.begin();
   while(it != groupid.end()) histogram[*it++]++;
   // Print out the frequencies of the values in v
   std::copy(histogram.begin(),histogram.end(),std::ostream_iterator< int >(std::cout, " "));
   std::cout << std::endl;
   // Find the mode
   int mode = std::max_element(histogram.begin(),histogram.end()) - histogram.begin();
   return mode;
   */
}



