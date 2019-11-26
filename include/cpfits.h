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

class CPFITS_BASE{
  
protected:
  CPFITS_BASE(){}

  CPFITS_BASE(CPFITS_BASE &&ff){
    fptr = ff.fptr;
    ff.fptr = nullptr;
  }
  
  void operator=(CPFITS_BASE &&ff){
    fptr = ff.fptr;
    ff.fptr = nullptr;
  }
  
  void check_status(int status,std::string s = "Error in CPFITS"){
    if(status){
      fits_report_error(stderr, status);
      std::cerr << s << std::endl;
      throw std::invalid_argument(s);
    }
  }
  fitsfile *fptr;
public:

  /////////////////////////////////
  // HDU-level Routines
  /////////////////////////////////
  
  /// returns the number of tabels or images
  int get_num_hdus(){
    int hdunum;
    int status = 0;
    fits_get_num_hdus(fptr, &hdunum, &status);
    check_status(status);
    return hdunum;
  }
  /// returns the current image / table number, starts with 1 not 0
  int get_current_hdu_num(){
    int hdunum;
    fits_get_hdu_num(fptr,&hdunum);
    return hdunum;
  }

  /// change the current table number
  int change_hdu(int i){
    int hdunum = get_num_hdus();
    if(i > hdunum || i < 1){
      throw std::invalid_argument("out of range");
    }
    int hdutype;
    int status = 0;
    fits_movabs_hdu(fptr,i,&hdutype,&status);
    check_status(status);
    reset_imageInfo();
 
    return hdutype;
  }
  
  int get_current_hdu_type(){
    int hdutype,status = 0;
    fits_get_hdu_type(fptr,&hdutype,&status);
    check_status(status);
    
    return hdutype;
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
  void reset_imageInfo(){
    
    int status = 0;
    fits_get_img_type(fptr,&bitpix,&status);
    check_status(status);
    fits_get_img_dim(fptr,&Ndims,&status);
    check_status(status);
    sizes.resize(Ndims);
    fits_get_img_size(fptr,Ndims,sizes.data(),&status);
    check_status(status);
    
    //std::cout << sizes.size() << std::endl;
  }
  void imageDimensions(std::vector<long> &my_sizes){
    my_sizes.resize(sizes.size());
    my_sizes = sizes;
  }

  /// number of keywords for current table, returns != 0 if keyname is not found
  int nkyewords(){
    int keysexist,status=0;
    return fits_get_hdrspace(fptr,&keysexist,NULL,&status);
  }
  
  /// read a key value for the current table / image, returns 0 if the key word does not exit
  int readKey(std::string keyname,double &value){
    int status = 0;
    return fits_read_key(fptr,TDOUBLE,keyname.c_str(),
                         &value,NULL,&status);
  }
  int readKey(std::string keyname,float &value){
    int status = 0;
    return fits_read_key(fptr,TFLOAT,keyname.c_str(),
                         &value,NULL,&status);
  }
  int readKey(std::string keyname,int &value){
    int status = 0;
    return fits_read_key(fptr,TINT,keyname.c_str(),
                         &value,NULL,&status);
  }
  int readKey(std::string keyname,size_t &value){
    int status = 0;
    return fits_read_key(fptr,TULONG,keyname.c_str(),
                         &value,NULL,&status);
  }
  int readKey(std::string keyname,long &value){
    int status = 0;
    return fits_read_key(fptr,TLONG,keyname.c_str(),
                         &value,NULL,&status);
  }

protected:
  int Ndims;
  std::vector<long> sizes;
  int bitpix;
};

/// read only fits file interface
class CPFITS_READ : public CPFITS_BASE {
private:
  // make it uncopyable
  //CPFITS_READ(CPFITS_READ &);
  //CPFITS_READ operator=(CPFITS_READ );
  
public:
  CPFITS_READ(std::string filename,bool verbose = false) {
    int status = 0;
    fits_open_file(&fptr,filename.c_str(), READONLY, &status);
    //    print any error messages
    check_status(status,"Problem with input fits file.");
    reset_imageInfo();

    
    if(verbose){
      std::cout << "Opening file : " << filename << std::endl;
      std::cout << "      status : " << status << std::endl;
    }
  }
  
  ~CPFITS_READ(){
    int status = 0;
    fits_close_file(fptr, &status);
    check_status(status,"Problem closing fits file!");
  }
  
  CPFITS_READ(CPFITS_READ &&ff):CPFITS_BASE(std::move(ff)){};
  
  void operator=(CPFITS_READ &&ff){
     CPFITS_BASE::operator=(std::move(ff));
  }
  
  /// close old and reopen a new file
  void reset(std::string filename){
    int status = 0;
    fits_close_file(fptr, &status);
    check_status(status);
    fits_open_file(&fptr,filename.c_str(), READONLY, &status);
    check_status(status);
    reset_imageInfo();
  }

  /// read the whole image into a vector
  int read(std::vector<double> &output,std::vector<long> &size){
    
    //int dtype;
    imageDimensions(size);
    long nelements=1;
    for(long n : size) nelements *= n;
    std::vector<long> start(size.size(),1);

    output.resize(nelements);
    int status = 0;
    return fits_read_pix(fptr,TDOUBLE,start.data(),nelements
                         ,NULL,output.data(),NULL, &status);
  }
  /// read the whole image into a vector
  int read(std::vector<float> &output,std::vector<long> &size){
    
    //int dtype;
    imageDimensions(size);
    long nelements=1;
    for(long n : size) nelements *= n;
    std::vector<long> start(size.size(),1);
    output.resize(nelements);
    
    int status = 0;
    return fits_read_pix(fptr,TFLOAT,start.data(),nelements
                         ,NULL,output.data(),NULL, &status);
  }
  /// read the whole image into a valarray
  int read(std::valarray<float> &output,std::vector<long> &size){
    
    //int dtype;
    imageDimensions(size);
    long nelements=1;
    for(long n : size) nelements *= n;
    std::vector<long> start(size.size(),1);

    output.resize(nelements);
    
    int status = 0;
    return fits_read_pix(fptr,TFLOAT,start.data(),nelements
                         ,NULL,&output[0],NULL, &status);
  }
  int read(std::valarray<double> &output,std::vector<long> &sizes){
    //std::cout << sizes.size() << std::endl;
    imageDimensions(sizes);
    std::cout << sizes.size() << std::endl;
    std::cout << sizes[0] << " " << sizes[1] << std::endl;
    long long nelements=1;
    for(long n : sizes) nelements *= n;
    std::vector<long> start(sizes.size(),1);
    output.resize(nelements);
    
    int status = 0;
    int error = fits_read_pix(fptr,TDOUBLE,start.data(),nelements
                            ,NULL,&output[0],NULL, &status);
    check_status(status,"PixMap Read Failure");

    return error;
  }

  /// read nelements in order from image to output array
  int read_block(double *output,long nelements,long * start){
    int status =0;
    return fits_read_pix(fptr,TDOUBLE,start,nelements
                         ,NULL,output,NULL, &status);
  }
  int read_block(float *output,long nelements,long * start){
    int status = 0;
    return fits_read_pix(fptr,TFLOAT,start,nelements
                         ,NULL,output,NULL, &status);
  }

  
  /// read a rectangular subset of the image
  int read_subset(double *output,long *lowerleft,long *upperright){
    long inc[2] = {1,1};  // this is a step
    int status = 0;
    return fits_read_subset(fptr,TDOUBLE,lowerleft,upperright,inc
                            ,NULL,output,NULL,&status);
  }
  /// read a rectangular subset of the image
  int read_subset(float *output,long *lowerleft,long *upperright){
    long inc[2] = {1,1};  // this is a step
    int status = 0;
    return fits_read_subset(fptr,TFLOAT,lowerleft,upperright,inc
                            ,NULL,output,NULL,&status);
  }

};

class CPFITS_WRITE : public CPFITS_BASE {
  
private:
  // make it uncopyable
  CPFITS_WRITE(CPFITS_READ &);
  CPFITS_WRITE operator=(CPFITS_WRITE );
  
public:
  CPFITS_WRITE(std::string filename,bool append = false,bool verbose = false) {
    
    if(!append){
      if(verbose) std::cout << "Creating file : " << filename << std::endl;
      if(filename[0] != '!') filename = "!" + filename;
      int status = 0;
      fits_create_file(&fptr, filename.c_str(), &status);
      check_status(status);
    }else{
      int status = 0;
      fits_open_file(&fptr,filename.c_str(),READWRITE,&status);
      reset_imageInfo();

      if(status==104){  // create file if it does not exist
        status = 0;
        fits_create_file(&fptr, filename.c_str(), &status);
        check_status(status);

        if(verbose) std::cout << "Creating file : " << filename << std::endl;
      }else{
        if(verbose) std::cout << "Appending to file : " << filename << std::endl;
        check_status(status);
      }
    }
  }
  
  ~CPFITS_WRITE(){
    int status = 0;
    fits_close_file(fptr, &status);
    check_status(status);
   }

  /// add a new image or cube to the file
  void write_image(std::vector<double> &im,std::vector<long> &size){
    size_t n=1;
    for(size_t i : size) n *= i;
    assert(im.size() == n);
    int status=0;
    fits_create_img(fptr,DOUBLE_IMG,size.size(),
                    size.data(), &status);
    check_status(status);
    std::vector<long> fpixel(size.size(),1);
    fits_write_pix(fptr,TDOUBLE,fpixel.data(),
                          im.size(),im.data(),&status);
    check_status(status);
  }
  void write_image(std::vector<float> &im,std::vector<long> &size){
    size_t n=1;
    for(size_t i : size) n *= i;
    assert(im.size() == n);
    int status = 0;
    fits_create_img(fptr,FLOAT_IMG,size.size(),
                    size.data(), &status);
    check_status(status);
    std::vector<long> fpixel(size.size(),1);
    fits_write_pix(fptr,TFLOAT,fpixel.data(),
                          im.size(),im.data(),&status);
    check_status(status);
  }
  void write_image(std::valarray<double> &im,std::vector<long> &size){
    size_t n=1;
    for(size_t i : size) n *= i;
    assert(im.size() == n);
    int status = 0;
    fits_create_img(fptr,DOUBLE_IMG,size.size(),
                    size.data(), &status);
    check_status(status);

    std::vector<long> fpixel(size.size(),1);
    fits_write_pix(fptr,TDOUBLE,fpixel.data(),
                          im.size(),&im[0],&status);
    check_status(status);
  }
  int write_image(std::valarray<float> &im,std::vector<long> &size){
    size_t n=1;
    for(size_t i : size) n *= i;
    assert(im.size() == n);
    int status = 0;
    fits_create_img(fptr,FLOAT_IMG,size.size(),
                    size.data(), &status);
    check_status(status);
    std::vector<long> fpixel(size.size(),1);
    fits_write_pix(fptr,TFLOAT,fpixel.data(),
                          im.size(),&im[0],&status);
    check_status(status);
  }

  /// add or replace a key value in the header of the current table / image
  int writeKey(std::string key,double value,std::string comment ){
    int status = 0;
    return fits_update_key(fptr,TDOUBLE,key.c_str(),
                           &value,comment.c_str(),&status);
  }
  int writeKey(std::string key,float value,std::string comment ){
    int status = 0;
    return fits_update_key(fptr,TFLOAT,key.c_str(),
                           &value,comment.c_str(),&status);
  }
  int writeKey(std::string key,int value,std::string comment ){
    int status = 0;
    return fits_update_key(fptr,TINT,key.c_str(),
                           &value,comment.c_str(),&status);
  }
  int writeKey(std::string key,long value,std::string comment ){
    int status = 0;
    return fits_update_key(fptr,TLONG,key.c_str(),
                           &value,comment.c_str(),&status);
  }
  int writeKey(std::string key,size_t value,std::string comment ){
    int status = 0;
    return fits_update_key(fptr,TULONG,key.c_str(),
                           &value,comment.c_str(),&status);
  }

};

#endif /* cpfits_h */
