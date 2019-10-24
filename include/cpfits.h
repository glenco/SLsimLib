//
//  cpfits.h
//  GLAMER
//
//  Created by Robert Benton Metcalf on 24/10/2019.
//
// This is a class that is meant to replace the CCFITS library
// so that it can be eliminated as a dependancy.
//
// These classes are for reading from and writing to fits files.
// Only
//
// for more infermation of cfitio routines see :
//   https://heasarc.gsfc.nasa.gov/docs/software/fitsio/quick/node9.html

#ifndef cpfits_h
#define cpfits_h

#include <string.h>
#include <stdio.h>
 #include "fitsio.h"

/// types of fits tables
enum HDUtype {IMAGE_HDU, ASCII_TBL, BINARY_TBL};

class COFITS_BASE{
  
protected:
  CPFITS_BASE()::status(0){}
  
  fitsfile *fptr;
  int status;
public:

  /////////////////////////////////
  // HDU-level Routines
  /////////////////////////////////
  
  /// returns the number of tabels
  int get_num_hdus(){
    int hdunum;
    return fits_get_num_hdus(fptr, &hdunum, &status);
  }
  /// returns the current table number, 1...
  int get_current_hdu_num(){
    int hdunum;
    return fits_get_hdu_num(fptr,&hdunum);
  }

  /// change the current table number
  HDUtype change_hdu(int i){
    int hdunum = get_num_hdus();
    if(i > hdum || i < 1){
      throw std::invalid
    }
    HDUtype hdutype;
    fits_movabs_hdu(fptr,i,&hdutype,&status);
    
    return hdutype;
  }
  
  HDUtype get_current_hdu_type(){
    HDUtype hdutype;
    return fits_get_hdu_type(fptr,&hdutype,&status)
  }
  
  /////////////////////////////////
  // image information of current HDU
  /////////////////////////////////

  /** \brief gives the data type and dimensions of the current table
   BYTE_IMG      =   8   ( 8-bit byte pixels, 0 - 255)
   SHORT_IMG     =  16   (16 bit integer pixels)
   LONG_IMG      =  32   (32-bit integer pixels)
   LONGLONG_IMG  =  64   (64-bit integer pixels)
   FLOAT_IMG     = -32   (32-bit floating point pixels)
   DOUBLE_IMG    = -64   (64-bit floating point pixels)
   */
  void imageInfo(
                 ,int &bitpix /// type of data
                 ,std::vector<long> &size
                 ){
    int ndim /// number of dimensions (keyword NAXIS)
    size.resize(10);
    fits_get_img_param(fptr,10,&bitpix,&ndim,size.data(),&status)
    size.resize(ndim);
  }

  /// number of keywords for current table, returns != 0 if keyname is not found
  int nkyewords(){
    int keysexist;
    return fits_get_hdrspace(fptr,&keysexist,NULL,&status);
  }
  int readKey(std::string &keyname,double &value){
    return fits_read_key(fptr,TDOUBLE,keyname.c_string,
                         &value,NULL,&status);
  }
  int readKey(std::string &keyname,float &value){
    return fits_read_key(fptr,TFLOAT,keyname.c_string,
                         &value,NULL,&status);
  }
  int readKey(std::string &keyname,int &value){
    return fits_read_key(fptr,TINT,keyname.c_string,
                         &value,NULL,&status);
  }
  int readKey(std::string &keyname,size_t &value){
    return fits_read_key(fptr,TULONG,keyname.c_string,
                         &value,NULL,&status);
  }

}

/// read only fits file interface
class CPFITS_READ : public CPFITS_BASE {
private:
  // make it uncopyable
  CPFITS_READ(CPFITS_READ &);
  CPFITS_READ operator=(CPFITS_READ );
  
public:
  CPFITS_READ(std::string filename,bool verboe = false) {
    fits_open_file(&fptr,filename.c, READONLY, &status);
    
    //    print any error messages
    if (status) fits_report_error(stderr, status);
    
    if(verbose){
      std::cout << "Opening file : " << filename << endl;
      std::cout << "      status : " << status << endl;
    }
  }
  
  ~CFITS_READ(){
    fits_close_file(fptr, &status);
    if (status) fits_report_error(stderr, status);
  }
  
  /// close old and reopen a new file
  void reset(std::string filename){
    fits_close_file(fptr, &status);
    assert(status==0);
    fits_open_file(&fptr,filename.c, READONLY, &status);
    //    print any error messages
    if (status) fits_report_error(stderr, status);
  }

  /// read the whole image into a vector
  int read(std::vector<double> &output){
    assert(start.size()==2);
    
    std::vector<long> size;
    int dtype;
    imageInfo(dtype,size);
    long nelements=1;
    for(long n : size) nelements *= n;
    std::vector<long> start(nelements,1);
    output.resize(nelements);
    
    return fits_read_pix(fptr,TDOUBLE,start.data(),nelements
                         ,NULL,output.data(),NULL, &status);
  }
  /// read the whole image into a vector
  int read(std::vector<float> &output){
    assert(start.size()==2);
    
    std::vector<long> size;
    int dtype;
    imageInfo(dtype,size);
    long nelements=1;
    for(long n : size) nelements *= n;
    std::vector<long> start(nelements,1);
    output.resize(nelements);
    
    return fits_read_pix(fptr,TFLOAT,start.data(),nelements
                         ,NULL,output.data(),NULL, &status);
  }

  /// read nelements in order from image to output array
  int read_block(double *output,long nelements,long * start){
    return fits_read_pix(fptr,TDOUBLE,start,nelements
                         ,NULL,output,NULL, &status);
  }
  int read_block(float *output,long nelements,long * start){
    return fits_read_pix(fptr,TFLOAT,start,nelements
                         ,NULL,output,NULL, &status);
  }

  
  /// read a rectangular subset of the image
  int read_subset(double *output,long *lowerleft,long *uppertright){
    long inc[2] = {1,1};  // this is a step
    return fits_read_subset(fptr,TDOUBLE,lowerleft,upperright,&inc
                            ,NULL,output,NULL,&status);
  }
  /// read a rectangular subset of the image
  int read_subset(float *output,long *lowerleft,long *uppertright){
    long inc[2] = {1,1};  // this is a step
    return fits_read_subset(fptr,TFLOAT,lowerleft,upperright,&inc
                            ,NULL,output,NULL,&status);
  }

}

class CPFITS_WRITE : public CPFITS_BASE {
  
private:
  // make it uncopyable
  CPFITS_WRITE(CPFITS_READ &);
  CPFITS_WRITE operator=(CPFITS_WRITE );
  
public:
  CPFITS_WRITE(std::string filename,bool verboe = false) {
    fits_create_file(&fptr, filename.c_string(), &status);
    
    
    //    print any error messages
    if (status) fits_report_error(stderr, status);
    
    if(verbose){
      std::cout << "Creating file : " << filename << endl;
      std::cout << "       status : " << status << endl;
    }
  }
  
  ~CFITS_WRITE(){
    fits_close_file(fptr, &status);
    if (status) fits_report_error(stderr, status);
  }

  /// add a new image or cube to the file
  int write_image(std::vector<double> &im,std::vector<long> &size){
    
    assert(im.size() == size[0]*siz[1]);
    fits_create_img(fptr,DOUBLE_IMG,size.size(),
                    size.data(), &status);
    std::vector<long> fpixel(size.sizs(),1);
    return fits_write_pix(fptr,TDOUBLE,fpixel.data(),
                          im.size(),im.data(),&status);
  }
  int writeImage(std::vector<float> &im,std::vector<long> &size){
    assert(im.size() == size[0]*siz[1]);
    fits_create_img(fptr,FLOAT_IMG,size.size(),
                    size.data(), &status);
    std::vector<long> fpixel(size.sizs(),1);
    return fits_write_pix(fptr,TFLOAT,fpixel.data(),
                          im.size(),im.data(),&status);
  }

  /// add or replace a key value in the header
  int writeKey(std::string &key,double value,std::string &comment ){
    return fits_update_key(fptr,TDOUBLE,key.c_string(),
                           &value,comment.c_string(),&status);
  }
  int writeKey(std::string &key,float value,std::string &comment ){
    return fits_update_key(fptr,TFLOAT,key.c_string(),
                           &value,comment.c_string(),&status);
  }
  int writeKey(std::string &key,int value,std::string &comment ){
    return fits_update_key(fptr,TINT,key.c_string(),
                           &value,comment.c_string(),&status);
  }
  int writeKey(std::string &key,long value,std::string &comment ){
    return fits_update_key(fptr,TLONG,key.c_string(),
                           &value,comment.c_string(),&status);
  }
  int writeKey(std::string &key,size_t value,std::string &comment ){
    return fits_update_key(fptr,TULONG,key.c_string(),
                           &value,comment.c_string(),&status);
  }

}

#endif /* cpfits_h */
