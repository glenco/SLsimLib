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
#include <iostream>
#include <stdio.h>
#include <valarray>
#include <map>
#include "fitsio.h"
#include <mutex>

#include <Tree.h>


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
  
  std::mutex mutex_lock;
public:

  /////////////////////////////////
  // HDU-level Routines
  /////////////////////////////////
  
  /// returns the number of tabels or images
  int get_num_hdus(){
    std::lock_guard<std::mutex> hold(mutex_lock);
    int hdunum;
    int status = 0;
    fits_get_num_hdus(fptr, &hdunum, &status);
    check_status(status);
    return hdunum;

  }
  /// returns the current image / table number, starts with 1 not 0
  int get_current_hdu_num(){
    std::lock_guard<std::mutex> hold(mutex_lock);
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
    {
      std::lock_guard<std::mutex> hold(mutex_lock);
      fits_movabs_hdu(fptr,i,&hdutype,&status);
    }
    check_status(status);
    reset_imageInfo();
 
    return hdutype;
  }
  
  int get_current_hdu_type(bool verbose=false){

    int hdutype,status = 0;
    fits_get_hdu_type(fptr,&hdutype,&status);
    check_status(status);
    
    if(verbose){
      std::cout << "type: " ;
      switch (hdutype) {
        case 0:
          std::cout << "Image" << std::endl;
          break;
        case 1:
          std::cout << "ASCII table" << std::endl;
          break;
        case 2:
          std::cout << "Binary table" << std::endl;
          break;

        default:
          break;
      }
    }
    
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
    std::lock_guard<std::mutex> hold(mutex_lock);

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
    std::lock_guard<std::mutex> hold(mutex_lock);
    int keysexist,status=0;
    return fits_get_hdrspace(fptr,&keysexist,NULL,&status);
  }
  
  /// read a key value for the current table / image, returns 0 if the key word does not exit
  int readKey(std::string keyname,double &value){
    std::lock_guard<std::mutex> hold(mutex_lock);
    int status = 0;
    return fits_read_key(fptr,TDOUBLE,keyname.c_str(),
                         &value,NULL,&status);
  }
  int readKey(std::string keyname,float &value){
    std::lock_guard<std::mutex> hold(mutex_lock);
    int status = 0;
    return fits_read_key(fptr,TFLOAT,keyname.c_str(),
                         &value,NULL,&status);
  }
  int readKey(std::string keyname,int &value){
    std::lock_guard<std::mutex> hold(mutex_lock);
    int status = 0;
    return fits_read_key(fptr,TINT,keyname.c_str(),
                         &value,NULL,&status);
  }
  int readKey(std::string keyname,size_t &value){
    std::lock_guard<std::mutex> hold(mutex_lock);
    int status = 0;
    return fits_read_key(fptr,TULONG,keyname.c_str(),
                         &value,NULL,&status);
  }
  int readKey(std::string keyname,long &value){
    std::lock_guard<std::mutex> hold(mutex_lock);
    int status = 0;
    return fits_read_key(fptr,TLONG,keyname.c_str(),
                         &value,NULL,&status);
  }
  int readHeader(std::vector<std::string> &headercards,bool verbose=false){
    std::lock_guard<std::mutex> hold(mutex_lock);
    int nkeys,status=0,hdupos;
    char card[FLEN_CARD];
    
    headercards.clear();
    fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */
    fits_get_hdrspace(fptr, &nkeys, NULL, &status); // get # of keywords
    
    if(verbose) printf("Header listing for HDU #%d:\n", hdupos);

    for (int ii = 1; ii <= nkeys; ii++) { // Read and print each keywords
      if (fits_read_record(fptr, ii, card, &status))break;
      if(verbose) printf("%s\n", card);
      headercards.push_back(card);
    }
    if(verbose) printf("END\n\n");  // terminate listing with END

      //fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
    return status;
  }

  /// print all header information
  int printHeader(){
    std::vector<std::string> cards;
    return readHeader(cards,true);
  }
  
   
  //int fits_parse_value(char *card, char *value, char *comment, int *status)
  //int fits_get_keytype(char *value, char *dtype, int *status)
  
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
  CPFITS_READ(std::string filename,std::string extension=""
              ,bool verbose = false) {

    int status = 0;
    //fits_open_file(&fptr,filename.c_str(), READONLY, &status);
    if(extension==""){
      fits_open_image(&fptr,filename.c_str(), READONLY, &status);
    }else{
      filename = filename + extension;
      fits_open_image(&fptr,filename.c_str(), READONLY, &status);
    }
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
    //check_status(status,"Problem closing fits file!");
  }
  
  CPFITS_READ(CPFITS_READ &&ff):CPFITS_BASE(std::move(ff)){};
  
  void operator=(CPFITS_READ &&ff){
     CPFITS_BASE::operator=(std::move(ff));
  }
  
  /// close old and reopen a new file
  void reset(std::string filename,std::string extension=""){
    std::lock_guard<std::mutex> hold(mutex_lock);
    int status = 0;
    fits_close_file(fptr, &status);
    check_status(status);
    //fits_open_image(&fptr,filename.c_str(), READONLY, &status);
    if(extension==""){
      fits_open_image(&fptr,filename.c_str(), READONLY, &status);
    }else{
      filename = filename + extension;
      fits_open_image(&fptr,filename.c_str(), READONLY, &status);
    }
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
    std::lock_guard<std::mutex> hold(mutex_lock);
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
    std::lock_guard<std::mutex> hold(mutex_lock);
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
    std::lock_guard<std::mutex> hold(mutex_lock);
    return fits_read_pix(fptr,TFLOAT,start.data(),nelements
                         ,NULL,&output[0],NULL, &status);
  }
  int read(std::valarray<double> &output,std::vector<long> &sizes){
    //std::cout << sizes.size() << std::endl;
    imageDimensions(sizes);
    //std::cout << sizes.size() << std::endl;
    //std::cout << sizes[0] << " " << sizes[1] << std::endl;
    long long nelements=1;
    for(long n : sizes) nelements *= n;
    std::vector<long> start(sizes.size(),1);
    output.resize(nelements);
    
    int status = 0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    int error = fits_read_pix(fptr,TDOUBLE,start.data(),nelements
                            ,NULL,&output[0],NULL, &status);
    check_status(status,"PixMap Read Failure");

    return error;
  }

  /// read the whole image into a vector
  template<typename T>
  int read(std::vector<std::vector<T> > &output){
  
    std::vector<long> size;
    std::vector<T> tmp;
    int status = read(tmp,size);
    output.resize(size[0]);
    size_t k = 0;
    for(size_t i=0 ; i < size[0] ; ++i){
      output[i].resize(size[1]);
      for(T &x : output[i]){
        x = tmp[k++];
      }
    }
    
    return status;
  }
  /// read the whole image into a vector
  
  /// read nelements in order from image to output array
  int read_block(double *output,long nelements,long * start){
    int status =0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    return fits_read_pix(fptr,TDOUBLE,start,nelements
                         ,NULL,output,NULL, &status);
  }
  int read_block(float *output,long nelements,long * start){
    int status = 0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    return fits_read_pix(fptr,TFLOAT,start,nelements
                         ,NULL,output,NULL, &status);
  }

  
  /// read a rectangular subset of the image
  int read_subset(double *output,long *lowerleft,long *upperright){
    long inc[2] = {1,1};  // this is a step
    int status = 0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    return fits_read_subset(fptr,TDOUBLE,lowerleft,upperright,inc
                            ,NULL,output,NULL,&status);
  }
  /// read a rectangular subset of the image
  int read_subset(float *output,long *lowerleft,long *upperright){
    long inc[2] = {1,1};  // this is a step
    int status = 0;
    std::lock_guard<std::mutex> hold(mutex_lock);
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
      fits_open_image(&fptr,filename.c_str(),READWRITE,&status);
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

   /// add a new array to the file
   void write_array(std::vector<double> &array){
     long n = array.size();
     int status=0;
     std::lock_guard<std::mutex> hold(mutex_lock);
     fits_create_img(fptr,DOUBLE_IMG,1,&n, &status);
     check_status(status);
     std::vector<long> fpixel(1,1);
     fits_write_pix(fptr,TDOUBLE,fpixel.data(),
                           n,array.data(),&status);
     check_status(status);
   }
   /// add a new array to the file
   void write_array(std::vector<float> &array){
     long n = array.size();
     int status=0;
     std::lock_guard<std::mutex> hold(mutex_lock);
     fits_create_img(fptr,FLOAT_IMG,1,&n, &status);
     check_status(status);
     std::vector<long> fpixel(1,1);
     fits_write_pix(fptr,TFLOAT,fpixel.data(),
                           n,array.data(),&status);
     check_status(status);
   }
  
  /// add a new image or cube to the file
  void write_image(std::vector<double> &im,std::vector<long> &size){
    size_t n=1;
    for(size_t i : size) n *= i;
    assert(im.size() == n);
    int status=0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    fits_create_img(fptr,DOUBLE_IMG,size.size(),
                    size.data(), &status);
    check_status(status);
    std::vector<long> fpixel(size.size(),1);
    fits_write_pix(fptr,TDOUBLE,fpixel.data(),
                          im.size(),im.data(),&status);
    check_status(status);
  }
  void write_image(double *im,std::vector<long> &size){
    size_t n=1;
    for(size_t i : size) n *= i;
    int status=0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    fits_create_img(fptr,DOUBLE_IMG,size.size(),
                    size.data(), &status);
    check_status(status);
    std::vector<long> fpixel(size.size(),1);
    fits_write_pix(fptr,TDOUBLE,fpixel.data(),
                          n,im,&status);
    check_status(status);
  }
 void write_image(std::vector<float> &im,std::vector<long> &size){
    size_t n=1;
    for(size_t i : size) n *= i;
    assert(im.size() == n);
    int status = 0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    fits_create_img(fptr,FLOAT_IMG,size.size(),
                    size.data(), &status);
    check_status(status);
    std::vector<long> fpixel(size.size(),1);
    fits_write_pix(fptr,TFLOAT,fpixel.data(),
                          im.size(),im.data(),&status);
    check_status(status);
  }
  void write_image(float *im,std::vector<long> &size){
    size_t n=1;
    for(size_t i : size) n *= i;
    int status=0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    fits_create_img(fptr,FLOAT_IMG,size.size(),
                    size.data(), &status);
    check_status(status);
    std::vector<long> fpixel(size.size(),1);
    fits_write_pix(fptr,TFLOAT,fpixel.data(),
                          n,im,&status);
    check_status(status);
  }
  void write_image(std::valarray<double> &im,std::vector<long> &size){
    size_t n=1;
    for(size_t i : size) n *= i;
    assert(im.size() == n);
    int status = 0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    fits_create_img(fptr,DOUBLE_IMG,size.size(),
                    size.data(), &status);
    check_status(status);

    status = 0;
    std::vector<long> fpixel(size.size(),1);
    fits_write_pix(fptr,TDOUBLE,fpixel.data(),
                          im.size(),&im[0],&status);
    check_status(status);
  }
  void write_image(std::valarray<float> &im,std::vector<long> &size){
    size_t n=1;
    for(size_t i : size) n *= i;
    assert(im.size() == n);
    int status = 0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    fits_create_img(fptr,FLOAT_IMG,size.size(),
                    size.data(), &status);
    check_status(status);
    std::vector<long> fpixel(size.size(),1);
    fits_write_pix(fptr,TFLOAT,fpixel.data(),
                          im.size(),&im[0],&status);
    check_status(status);
  }

  /// write a two dimensional image that is stored in a vector of vectors
  template<typename T>
  void write_image(std::vector<std::vector<T> > &im){
    std::vector<long> size(2);
    size[0] = im.size();
    size[1] = im[0].size();
    
    size_t k=0;
    std::vector<T> tmp(size[0]*size[1]);
    for(long i=0; i<size[0] ; ++i){
      for(long j=0; j<size[1] ; ++j){
        tmp[k++] = im[i][j];
      }
    }
    
    write_image(tmp,size);
  }

  /// add or replace a key value in the header of the current table / image
  int writeKey(std::string key,std::string value,std::string comment ){
    int status = 0;
    char *s = strdup(value.c_str());
    std::lock_guard<std::mutex> hold(mutex_lock);
    int error = fits_write_key(fptr,TSTRING,key.c_str(),
                           s,comment.c_str(),&status);
    delete[] s;
    return error;
  }
  int writeKey(std::string key,double value,std::string comment ){
    int status = 0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    return fits_update_key(fptr,TDOUBLE,key.c_str(),
                           &value,comment.c_str(),&status);
  }
  int writeKey(std::string key,float value,std::string comment ){
    int status = 0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    return fits_update_key(fptr,TFLOAT,key.c_str(),
                           &value,comment.c_str(),&status);
  }
  int writeKey(std::string key,int value,std::string comment ){
    int status = 0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    return fits_update_key(fptr,TINT,key.c_str(),
                           &value,comment.c_str(),&status);
  }
  int writeKey(std::string key,long value,std::string comment ){
    int status = 0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    return fits_update_key(fptr,TLONG,key.c_str(),
                           &value,comment.c_str(),&status);
  }
  int writeKey(std::string key,size_t value,std::string comment ){
    int status = 0;
    std::lock_guard<std::mutex> hold(mutex_lock);
    return fits_update_key(fptr,TULONG,key.c_str(),
                           &value,comment.c_str(),&status);
  }

  /// write many header cards at once
  int writeHeader(std::vector<std::string> &cards){
    int status = 0;
    for(auto &str : cards){
      fits_write_record(fptr, str.c_str(), &status);
    }
    return status;
  }
  
  /// copy all the header information from the input file to the output file
  void copyHeader(CPFITS_BASE &input){
    std::vector<std::string> cards;
    input.readHeader(cards);
    writeHeader(cards);
  }

};

/// read only fits file interface
class CPFITS_READ_TABLES : public CPFITS_BASE {
private:
  // make it uncopyable
  //CPFITS_READ(CPFITS_READ &);
  //CPFITS_READ operator=(CPFITS_READ );
  
  long Nrow;  // current number of rows
  int Ncol;  // current number of columns
  
  void reset_tableInfo(){
    int status = 0;
    fits_get_num_rows(fptr, &Nrow, &status);
    fits_get_num_cols(fptr, &Ncol, &status);
  }
  std::vector<std::string> col_names;
  
  int get_colnames(std::vector<std::string> &names){
    
    int status=0;
    char colname[100];
    std::string c = "*";
    names.resize(Ncol);
    int colnum;
    
    //for(int i=0;i<Ncol;++i){
    //for(int i=0;i<10;++i){
    for(auto &name : names){
      fits_get_colname(fptr,CASESEN,(char*)(c.c_str()),colname, &colnum, &status);
      //names[i] = colname;
      name = colname;
      //std::cout << i << " " << names[i] << "  " << status << std::endl;
    }
    
    return status;
  }
public:
  
  CPFITS_READ_TABLES(std::string filename,int hdunum=2,bool verbose = false) {

    int status = 0;
    fits_open_table(&fptr,filename.c_str(), READONLY, &status);
    check_status(status,"Problem with input fits file.");
    if(hdunum>0){
      int num;
      fits_get_num_hdus(fptr, &num, &status);
      if(hdunum > num){
        throw std::runtime_error("Not enough HDUs in " + filename);
      }
      int hdutype;
      fits_movabs_hdu(fptr,hdunum, &hdutype,  &status);
      check_status(status,"Problem with input fits file. HDU "+std::to_string(hdunum));
      if(hdutype != BINARY_TBL){
        std::cerr << "HDU " << hdunum << " in " << filename << " is not a binary table" << std::endl;
        throw std::runtime_error("wrong HDUs in " + filename);
      }
    }

    reset_tableInfo();
    get_colnames(col_names);
    
    if(Nrow <= 0){
      std::cerr << "No data was found in file " << filename << std::endl;
    }
    
    if(verbose){
      std::cout << "Opening file : " << filename << std::endl;
      std::cout << Ncol << " columns" << std::endl;
      std::cout << Nrow << " rows" << std::endl;
    }
  }
  
  ~CPFITS_READ_TABLES(){
    int status = 0;
    fits_close_file(fptr, &status);
  }
  
  CPFITS_READ_TABLES(CPFITS_READ &&ff):CPFITS_BASE(std::move(ff)){};
  
  void operator=(CPFITS_READ_TABLES &&ff){
     CPFITS_BASE::operator=(std::move(ff));
  }
  
  long rows(){return Nrow;}
  long ncol(){return Ncol;}
  
  int column_index(std::string &colname){
    int col, status=0;
    
    fits_get_colnum(fptr, CASEINSEN,(char*)(colname.c_str())
    , &col, &status);
    check_status(status,"Column does not exist.");
    
    return col;
  }

  std::vector<std::string> get_colnames(){return col_names;}
  
  int read_column(std::string &colname,std::vector<float> &vecf){
   
    int col, status=0,anyval=0;
    col = column_index(colname);

    if(vecf.size() < Nrow ) vecf.resize(Nrow);
    
    fits_read_col_flt(fptr, col, 1, 1,Nrow,0,
                      vecf.data(),&anyval, &status);

    return status;
  }
  int read_column(std::string &colname,long firstsrow,long nelements,std::vector<float> &vecf,float nulval = -100){
   
    int status=0,anyval=0;
    int colnum = column_index(colname);
    
    if(vecf.size() < nelements ) vecf.resize(nelements);
    fits_read_col(fptr,TFLOAT,colnum, firstsrow, 1,nelements,&nulval,
                      vecf.data(),&anyval, &status);
    check_status(status);

    return status;
  }
  int read_column(int colnum,long firstsrow,long nelements,std::vector<float> &vecf,float nulval = -100){
   
    int status=0,anyval=0;
    
    if(vecf.size() < nelements ) vecf.resize(nelements);
    
    fits_read_col(fptr,TFLOAT,colnum, firstsrow, 1,nelements,&nulval,
                      vecf.data(), &anyval, &status);

    check_status(status);

    return status;
  }
  
  int read_column(std::string &colname,long firstsrow,long nelements,std::vector<double> &vecf,double nulval = -100){
   
    int status=0,anyval=0;
    int colnum = column_index(colname);
    
    if(vecf.size() < nelements ) vecf.resize(nelements);
    fits_read_col(fptr,TDOUBLE,colnum, firstsrow, 1,nelements,&nulval,
                      vecf.data(),&anyval, &status);
    check_status(status);

    return status;
  }
  int read_column(int colnum,long firstsrow,long nelements,std::vector<double> &vecf,double nulval = -100){
   
    int status=0,anyval=0;
    
    if(vecf.size() < nelements ) vecf.resize(nelements);
    fits_read_col(fptr,TDOUBLE,colnum, firstsrow, 1,nelements,&nulval,
                      vecf.data(),&anyval, &status);
    check_status(status);

    return status;
  }
  int read_column(std::string &colname,long firstsrow,long nelements,std::vector<long> &vecf,int nulval = -100){
   
    int status=0,anyval=0;
    int colnum = column_index(colname);
    
    if(vecf.size() < nelements ) vecf.resize(nelements);
    fits_read_col(fptr,TLONG,colnum, firstsrow, 1,nelements,&nulval,
                      vecf.data(),&anyval, &status);
    check_status(status);

    return status;
  }
  int read_column(int colnum,long firstsrow,long nelements,std::vector<long> &vecf,int nulval = -100){
   
    int status=0,anyval=0;
    
    if(vecf.size() < nelements ) vecf.resize(nelements);
    fits_read_col(fptr,TLONG,colnum, firstsrow, 1,nelements,&nulval,
                      vecf.data(),&anyval, &status);
    check_status(status);

    return status;
  }


  /// get data type of column
  int gettype(int colnum){
    int typecode,status=0;
    long repeat,width;
    fits_get_coltype(fptr,colnum,&typecode,
                     &repeat,&width,&status);
    check_status(status);
    
    return typecode;
  }
  
  int read_column(int column_index,std::vector<float> &vecf){
   
    int status=0,anyval=0;

    if(column_index < 1){
      std::cerr << "CPFITS: column number must be > 0" << std::endl;
      throw std::invalid_argument("Error in CPFITS");
    }
    if(column_index > Ncol){
      std::cerr << "CPFITS: there are only " << Ncol << " columns." << std::endl;
      throw std::invalid_argument("Error in CPFITS");
    }
    if(vecf.size() < Nrow ) vecf.resize(Nrow);
    
    fits_read_col_flt(fptr, column_index, 1, 1,Nrow,0,
                      vecf.data(),&anyval, &status);

    return status;
  }

  
  template <typename T>
  void read(std::vector<T> data,int ncol){
    
      int status = 0,anyval=0;
    long i, n, nread, ntodo, nrest;
    long buffer_size = 20000;
    std::vector<float> fitscol(buffer_size);
    
    for(int col = 1 ; col != ncol ; ++col){
      nread = 0;
      nrest = n;
      while (nrest) {
        //ntodo = MIN(nrest,buffer_size);
        ntodo = (nrest < buffer_size) ? nrest : buffer_size;
       fits_read_col_flt(fptr, col, nread + 1, 1, ntodo, 0,
                             fitscol.data(),&anyval, &status);
       for (i = 0; i < ntodo; i++) data[i + nread][col-1] = fitscol[i];
       nread += ntodo;
       nrest -= ntodo;
      }
    }
  }
};

/*** \brief Data frame for reading fits tables and minipulating the results.
 
 Missing entries are replaced with -100
 */

  template< typename T>
  class DataFrameFits{
 
  private:
    std::map<std::string,int> datamap;
    std::vector<std::string> all_column_names;
    std::vector<std::string> used_column_names;
    std::string filename;
    CPFITS_READ_TABLES cpfits;
    std::vector<std::vector<T> > data;
    long n0;
    std::vector<int> column_index;
  
  public:
    /// The data frame is attached to a file, but the constructor will not read all the data in the table.  Reading needs to be done with a seporate function.
    DataFrameFits(
                  std::string datafile   /// input catalog file in fits format
                  ,std::vector<std::string> &columns  /// if empty all columns are read and this will contain thier names, if not, only the listed columns are read
                  ,int hdu = 2
                  ,bool verbose = false
    ):filename(datafile),cpfits(filename,hdu,verbose),n0(1){
   
      all_column_names = cpfits.get_colnames();
      
      if(columns.size() == 0){
        used_column_names = all_column_names;
      }else{
        used_column_names = columns;
      }
      
      std::cout << "Reading " << used_column_names.size() << " columns." << std::endl;
      
      for(std::string a : used_column_names){
        column_index.push_back(cpfits.column_index(a));
      }
      
      for(int i=0 ; i<used_column_names.size() ; ++i){
        datamap[used_column_names[i]] = i;
      }
 
      data.resize(column_index.size());
    };
   
    /// returns all the column names in the fits table
    std::vector<std::string> get_all_columnnames(){return all_column_names;}
    std::vector<std::string> get_used_columnnames(){return all_column_names;}
    
    // clear data
    void clear(){
      n0=1;
      int ncol=data.size();
      for(int i=0; i<ncol ; ++i){
        data[i].clear();
      }
    }
 
    /** \brief reads in only the rows that meet the criterion given by accept criterion.
        accept -  a vector of functions that define acceptence in the order of the column names provided to the constructor
        add=true  - add more rows with a new acceptence criterion searched from the top of the file
        add=false - clear contents and start from where it last left off in the file
        maxsize -  maximum number of rows added, can read in batches by running again with add=false, default is the total number of rows in the file
     */
    int read(
             std::vector<std::function<bool(T &)> > &accept
             ,bool add=false
             ,long maxsize = -1
             ){
      long chunksize = 10000;
      int ncol = column_index.size();
      long nrow = cpfits.rows();
      int requirements = accept.size();

      if(add){
        n0=1;
      }else{
        for(int i=0; i<ncol ; ++i){
          data[i].clear();
        }
      }
  
      if(maxsize > 0) nrow = MIN(maxsize + n0,nrow);
      
      if(requirements > ncol){
        throw std::invalid_argument("Too many requirments");
      }
  
      std::vector<std::vector<T> > tdata(ncol);
      while(n0 < nrow){
        chunksize = MIN(nrow-n0+1,chunksize);
        
        if(data[0].capacity() - data[0].size() < chunksize ){
          for(int i=0; i<ncol ; ++i){
            data[i].reserve(chunksize + data[i].capacity());
          }
        }
      
        for(int i=0 ; i< ncol ; ++i){
          cpfits.read_column(column_index[i],n0,chunksize,tdata[i]);
        }
        n0 += chunksize;
      
        for(long j = 0 ; j < chunksize ; ++j){
          bool accpt = true;
          for(int i=0; i<requirements; ++i){
            accpt *= accept[i](tdata[i][j]);
          }
          
          if(accpt){
            for(int i=0; i<ncol ; ++i){
              data[i].push_back( tdata[i][j] ) ;
            }
          }
  
        }
      }
      
      return 1;
    };
    
    size_t read(
             std::function<bool(T,T)> &binary_accept
             ,std::pair<int,int> index_binary                  /// the columns to be compared
             ,std::vector<std::function<bool(T &)> > &unary_accept
             ,std::vector<int> index_unary
             ,bool add=false
             ,long maxsize = -1
             ){
      long chunksize = 10000;
      int ncol = column_index.size();
      long nrow = cpfits.rows();
      
      if(unary_accept.size() != index_unary.size() ) throw std::invalid_argument("index_unary wrong sized ");
 
      if(add){
        n0=1;
      }else{
        for(int i=0; i<ncol ; ++i){
          data[i].clear();
        }
      }
  
      if(maxsize > 0) nrow = MIN(maxsize + n0,nrow);
      
      if(index_unary.size() > ncol){
        throw std::invalid_argument("Too many requirments");
      }
  
      std::vector<std::vector<T> > tdata(ncol);
      while(n0 < nrow){
        chunksize = MIN(nrow-n0+1,chunksize);
        
        if(data[0].capacity() - data[0].size() < chunksize ){
          for(int i=0; i<ncol ; ++i){
            data[i].reserve(chunksize + data[i].capacity());
          }
        }
      
        for(int i=0 ; i< ncol ; ++i){
          cpfits.read_column(column_index[i],n0,chunksize,tdata[i]);
        }
        n0 += chunksize;
      
        for(long j = 0 ; j < chunksize ; ++j){
    
          bool accpt = binary_accept(tdata[index_binary.first][j],tdata[index_binary.second][j] );
          
          int k=0;
          for(int i : index_unary){
            accpt *= unary_accept[k++]( tdata[i][j] );
          }
          
          if( accpt ){
            for(int i=0; i<ncol ; ++i){
              data[i].push_back( tdata[i][j] ) ;
            }
          }
  
        }
      }
      
      return data.size();
    }

    /// read the next maxsize rows.  This will errase the rows already read
    int read(long maxsize=-1){
      int ncol = column_index.size();
      long nrow = cpfits.rows();
      if(maxsize==-1) maxsize = nrow;
      
      long chunksize = nrow-n0+1;
      //chunksize = MIN(maxsize,chunksize);
      chunksize = (maxsize < chunksize) ? maxsize : chunksize;
      
      for(int i=0 ; i< ncol ; ++i){
        cpfits.read_column(column_index[i],n0,chunksize,data[i]);
      }
      n0 += chunksize;
      
      return 1;
    }

    /// Find the ranges for each used column
    std::vector<std::vector<T> > ranges(bool verbose){
      std::vector<std::vector<T> > ranges;
      int ncol = column_index.size();
      long nrow = cpfits.rows();
      long chunksize = 100000;
      
      ranges.resize(ncol);
      for(auto &v : ranges) v.resize(2);
    
      std::vector<T> tdata;
      for(int i=0 ; i< ncol ; ++i){
        cpfits.read_column(column_index[i],1,1,tdata);
        for(T a : tdata){
          ranges[i][0] = a;
          ranges[i][1] = a;
        }
      }

      
      long n=2;
      while(n < nrow){

        //chunksize = MIN(nrow-n+1,chunksize);
        chunksize = (nrow-n+1 < chunksize) ? nrow-n+1 : chunksize;

      
        for(int i=0 ; i< ncol ; ++i){
          cpfits.read_column(column_index[i],n,chunksize,tdata);
          for(T a : tdata){
            ranges[i][0] = MIN(a,ranges[i][0]);
            ranges[i][1] = MAX(a,ranges[i][1]);
          }
        }
        n += chunksize;
      }
      
      if(verbose){
        for(int i=0 ; i< ncol ; ++i){
          std::cout << ranges[i][0] << " <= " << used_column_names[i] << " <= " << ranges[i][1] << std::endl;
        }
      }
      
      return ranges;
    }
    
    /// returns column by name
    std::vector<T>& operator[](const std::string &label){
      if(datamap.find(label) == datamap.end()){
        std::cerr << "No label - " << label << " - in " << filename <<std::endl;
        throw std::invalid_argument("no label");
      }
      return data[datamap[label]];
      for(auto c : all_column_names ) std::cout << c << " ";
      std::cout << std::endl;
      throw std::invalid_argument(label + " was not one of the columns of the galaxy data file :" + filename);
    };
    
    /// returns column by number in order of input column names or the original table in not provided
    std::vector<T>& operator[](int i){
      return data[i];
    };
    
    // sort by one of the columns in ascending order
    void sortby(std::string name){
      std::vector<size_t> index(data[0].size());
      size_t N = index.size();
      for(size_t i=0 ; i<N ; ++i) index[i] = i;
      
      Utilities::sort_indexes(data[datamap[name]],index);
      std::vector<T> tmp_v(N);
      
      for(size_t j=0 ; j<data.size() ; ++j){
        for(size_t i=0 ; i<N ; ++i){
          tmp_v[i] = data[j][index[i]];
        }
        swap(data[j],tmp_v);
      }
    }
    
    /// add a new column of zeros
    void addColumn(std::string name){
      int nc =  number_of_columns();
      long nr = number_of_rows();
      data.emplace_back(nr);
      all_column_names.push_back(name);
      used_column_names.push_back(name);
      datamap[name] = nc;
      //column_index.push_back(nc);
    }
    
    size_t number_of_rows(){return data[0].size();}
    size_t number_of_columns(){return data.size();}

};
#endif /* cpfits_h */
