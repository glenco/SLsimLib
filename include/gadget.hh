#ifndef GADGET_HPP
#define GADGET_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <sys/stat.h>
#include <sstream>
#include <string>
#include <cstring>
#include <math.h>
#include <vector>

#include "particle_types.h"

#define int4bytes int

using namespace std;

template <typename T>
constexpr std::size_t size_of(T const&) {
  return sizeof(T);
}
//class gadget: public gensim {
template<typename PType>
class GadgetFile {
  /* this is a specialized version of gensim which deals with GADGET files */
public:
  
  bool multipleFiles;
  int numfiles;
  size_t ntot=0;
  double time,redshift,omega,lambda,hubble;
  std::string filebasename;
  
  GadgetFile (string inpfn,std::vector<PType> &data);
  ~GadgetFile(){};
  
  void checkMultiple ();
  void openFile ();
  void closeFile ();
  
  void readBlock (const char *blockname);
  
  int find_block(FILE *fd,const char *label);
  
  int read_gadget_head(int *npart,double *massarr,double *time,double *redshift,double *omega, double *lambda, double *hubble, FILE *fd);
  
  bool FileExists(string strFilename);

  std::vector<PType> &p_data;

  int npart[6];
  double masstab[6];

private:
  
  struct pvector{
    float x,y,z;
  };
  
  ifstream::pos_type size;
  int4bytes blksize, swap;

  //float *rho, *mass, *u, *den, *pot;
  //int *ptype;
  //pvector *pos,*b,*vel;
  FILE *fd;
  int filecnt; 
  int np_file_start, np_file_end;
  
  int read_gadget_float(float *data,const char *label,FILE *fd);
  
  int read_gadget_float3(float *data,const char *label,FILE *fd);
  size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
  void swap_Nbyte(char *data,int n,int m);
};

template<typename PType>
GadgetFile<PType>::GadgetFile(string inpfn,std::vector<PType> &data):
multipleFiles(false),numfiles(0),filebasename(inpfn),p_data(data),swap(0)
,filecnt(0),np_file_start(0),np_file_end(-1)
{
  checkMultiple();
}

template<typename PType>
void GadgetFile<PType>::checkMultiple(){
  
  char *filename = new char [filebasename.size()+1];
  int n;
  filename[filebasename.size()]=0;
  memcpy(filename,filebasename.c_str(),filebasename.size());
  
  bool exists=FileExists(filename);
  if (exists) {cout << "The file exists..." << endl;
    
    if(!(fd = fopen(filename,"r")))
    {
      printf("Cant open file <%s> !\n",filename);
      exit(2);
    }else{
      cout << "Reading header" << endl;
      
      /*----------- READ HEADER TO GET GLOBAL PROPERTIES -------------*/
      n = read_gadget_head(npart,masstab,&time,&redshift,&omega,&lambda,&hubble,fd);
      
      ntot=0;
      for(int i=0;i<6;i++)
      {
        printf("PartSpecies %d, anz=%d, masstab=%f\n",i,npart[i],masstab[i]);
        ntot += npart[i];
      }
      printf("Omage %f Lambda %f hubble %f\n",omega,lambda,hubble);
      //numfiles = 1;
    }
    fclose(fd);
  } else {
    cout << "The file does not exist or the simulation is splitted in multiple files..." << endl;
    cout << "Trying multiple files..." << endl;
    int nmaxf=10;
    
    for (int i=0; i<=nmaxf-1; i++) {
      
      std::ostringstream filenum;
      string odot=".";
      filenum << i;
      string filebasename_=filebasename+odot+filenum.str();
      char *filename_ = new char [filebasename.size()+1];
      filename_[filebasename_.size()]=0;
      memcpy(filename_,filebasename_.c_str(),filebasename_.size());
      bool exists=FileExists(filename_);
      if (exists) {numfiles++;}
      
    }
    if (numfiles>0) {
      multipleFiles=true;
      cout << "The simulation is splitted into " << numfiles << " files" << endl;
      ntot=0;
      for (int j=0; j<=numfiles-1; j++) {
        std::ostringstream filenum;
        string odot=".";
        filenum << j;
        string filebasename_=filebasename+odot+filenum.str();
        char *filename_ = new char [filebasename.size()+1];
        filename_[filebasename_.size()]=0;
        memcpy(filename_,filebasename_.c_str(),filebasename_.size());
        
        if(!(fd = fopen(filename_,"r")))
        {
          printf("Cant open file <%s> !\n",filename_);
          exit(2);
        }
        else
        {
          cout << "-> Reading header# " << j << endl;
          
          /*----------- READ HEADER TO GET GLOBAL PROPERTIES -------------*/
          n = read_gadget_head(npart,masstab,&time,&redshift,&omega,&lambda,&hubble,fd);
          
          for(int i=0;i<6;i++)
          {
            printf("* PartSpecies %d, anz=%d, masstab=%f\n",i,npart[i],masstab[i]);
            ntot += npart[i];
          }
        }
      }
      fclose(fd);
    }
  }
  cout << "TOTAL NUMBER OF PARTICLES: " << ntot << endl;
  //now create a vector of particles
  printf("Allocating memory for a total of %zu particles...\n",ntot);
 
  p_data.resize(ntot);
  
  // pv.reserve(ntot);
  // for (int k=0; k<ntot; k++){
  //   particle* temp= new particle();
  //   pv.push_back(temp);
  // }
  
  cout << "done" << endl;
}

//*******************************************************************************
//*******************************************************************************
//*******************************************************************************
//*******************************************************************************
//*******************************************************************************
//*******************************************************************************


#define SKIP  {my_fread(&blksize,sizeof(int),1,fd); swap_Nbyte((char*)&blksize,1,4);}

template<typename PType>
void GadgetFile<PType>::openFile () {
  
  int n;
  int nstart;
  
  fd=0;
  
  if (numfiles==0) {
    
    char *filename = new char [filebasename.size()+1];
    filename[filebasename.size()]=0;
    memcpy(filename,filebasename.c_str(),filebasename.size());
    fd = fopen(filename,"r");
    
    /*----------- RED HEADER TO GET GLOBAL PROPERTIES -------------*/
    n = read_gadget_head(npart,masstab,&time,&redshift,&omega,&lambda,&hubble,fd);
    
    
    // now setting the particle type
    nstart=0;
    np_file_start=nstart;
    
    for (int ips=0; ips<6; ips++){
      if (npart[ips]>0){
        for (int i=nstart; i<nstart+npart[ips]; i++){
          p_data[i].type = ips;
        }
        nstart=nstart+npart[ips];
      }
    }
    np_file_end=np_file_start+ntot-1;
  } else {
    
    // if there are multiple files, these have to be opened one-by-one. Particles must be handled carefully: must define
    // np_file_* for each file and must update the number of the opened file
    
    // step 1: construct the file name. This is determined by the basename + the file number (which is updated afterwords)
    // NB: EVERY TIME THE FUNCTION IS CALLED A NEW FILE IS OPENED!
    
    
    std::ostringstream filenum;
    string odot=".";
    filenum << filecnt;
    string filebasename_=filebasename+odot+filenum.str();
    char *filename_ = new char [filebasename.size()+1];
    filename_[filebasename_.size()]=0;
    memcpy(filename_,filebasename_.c_str(),filebasename_.size());
    
    fd = fopen(filename_,"r");
    
    /*----------- RED HEADER TO GET GLOBAL PROPERTIES -------------*/
    n = read_gadget_head(npart,masstab,&time,&redshift,&omega,&lambda,&hubble,fd);
    int np_in_file;
    np_in_file=0;
    // now setting the particle type
    
    np_file_start=np_file_end+1;
    nstart=np_file_start;
    
    for (int ips=0; ips<6; ips++){
      if (npart[ips]>0){
        for (int i=nstart; i<nstart+npart[ips]; i++){
          //pv[i]->setPtype(ips);
          p_data[i].type = ips;
        }
        nstart=nstart+npart[ips];
        np_in_file=np_in_file+npart[ips];
      }
    }
    np_file_end=np_file_start+np_in_file-1;
    filecnt++;
  }
  //cout << "PARTICLES FROM " << np_file_start << " TO " << np_file_end << endl;
}

template<typename PType>
void GadgetFile<PType>::closeFile () {
  fclose(fd);
}

/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/
/*---------------------------- High Level Routines ----------------------------*/
/*-----------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to read the header information ---------------*/
/*-------- int *npart:    List of Particle numbers for spezies [0..5] ---------*/
/*-------- int *massarr:  List of masses for spezies [0..5] -------------------*/
/*-------- int *time:     Time of snapshot ------------------------------------*/
/*-------- int *massarr:  Redshift of snapshot --------------------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns number of read bytes ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
template<typename PType>
int GadgetFile<PType>::read_gadget_head(int *npart,double *massarr,double *time,double *redshift, double *omega, double *lambda, double *hubble, FILE *fd)
{
  int blocksize,dummysize;
  
  int dummyint;
  int dummyarr[6];
  double boxl;
  
  blocksize = find_block(fd,"HEAD");
  
  
  if(blocksize <= 0)
  {
    printf("Block <%s> not fond !\n","HEAD");
    exit(5);
  }
  else
  {
    dummysize=blocksize - 6 * sizeof(int) - 8 * sizeof(double);
    std::cout << blocksize << " " << dummysize << std::endl;
    SKIP;
    my_fread(npart,6*sizeof(int), 1, fd);        swap_Nbyte((char*)npart,6,4);
    my_fread(massarr,6*sizeof(double), 1, fd);   swap_Nbyte((char*)massarr,6,8);
    my_fread(time,sizeof(double), 1, fd);        swap_Nbyte((char*)time,1,8);
    my_fread(redshift,sizeof(double), 1, fd);    swap_Nbyte((char*)redshift,1,8);
    my_fread(&dummyint,sizeof(int), 1, fd);       swap_Nbyte((char*)&dummyint,1,4);
    my_fread(&dummyint,sizeof(int), 1, fd);       swap_Nbyte((char*)&dummyint,1,4);
    my_fread(dummyarr,6*sizeof(int), 1, fd);     swap_Nbyte((char*)dummyarr,6,4);
    my_fread(&dummyint,sizeof(int), 1, fd);       swap_Nbyte((char*)&dummyint,1,4);
    my_fread(&dummyint,sizeof(int), 1, fd);       swap_Nbyte((char*)&dummyint,1,4);
    my_fread(&boxl,sizeof(double), 1, fd);        swap_Nbyte((char*)&boxl,1,8);
    //my_fread(&omega_tmp,sizeof(double), 1, fd);   swap_Nbyte((char*)&omega_tmp,1,8);
    //my_fread(&lambda_tmp,sizeof(double), 1, fd);  swap_Nbyte((char*)&lambda_tmp,1,8);
    //my_fread(&hubble_tmp,sizeof(double), 1, fd);  swap_Nbyte((char*)&hubble_tmp,1,8);
    my_fread(omega,sizeof(double), 1, fd);   swap_Nbyte((char*)omega,1,8);
    my_fread(lambda,sizeof(double), 1, fd);  swap_Nbyte((char*)lambda,1,8);
    my_fread(hubble,sizeof(double), 1, fd);  swap_Nbyte((char*)hubble,1,8);
    fseek(fd,dummysize,1);
    SKIP;
  }
  return(blocksize);
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to swap ENDIAN -------------------------------*/
/*-------- char *data:    Pointer to the data ---------------------------------*/
/*-------- int n:         Number of elements to swap --------------------------*/
/*-------- int m:         Size of single element to swap ----------------------*/
/*--------                int,float = 4 ---------------------------------------*/
/*--------                double    = 8 ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
template<typename PType>
void GadgetFile<PType>::swap_Nbyte(char *data,int n,int m)
{
  int i,j;
  char old_data[16];
  
  if(swap>0)
  {
    for(j=0;j<n;j++)
    {
      memcpy(&old_data[0],&data[j*m],m);
      for(i=0;i<m;i++)
      {
        data[j*m+i]=old_data[m-i-1];
      }
    }
  }
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine find a block in a snapfile -------------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- char *label:   4 byte identifyer for block -------------------------*/
/*-------- returns length of block found, -------------------------------------*/
/*-------- the file fd points to starting point of block ----------------------*/
/*-----------------------------------------------------------------------------*/
template<typename PType>
int GadgetFile<PType>::find_block(FILE *fd,const char *label)
{
  int4bytes blocksize=0;
  char blocklabel[5]={"    "};
  
  rewind(fd);
  
  while(!feof(fd) && blocksize == 0)
  {
    SKIP;
    if(blksize == 134217728)
    {
#ifdef MY_DEBUG
      printf("Enable ENDIAN swapping !\n");
#endif
      swap=1-swap;
      swap_Nbyte((char*)&blksize,1,4);
    }
    if(blksize != 8)
    {
      printf("incorrect format (blksize=%d)!\n",blksize);
      exit(1);
    }
    else
    {
      my_fread(blocklabel, 4*sizeof(char), 1, fd);
      my_fread(&blocksize, sizeof(int4bytes), 1, fd);
      swap_Nbyte((char*)&blocksize,1,4);
#ifdef MY_DEBUG
      printf("Found Block <%s> with %d bytes\n",blocklabel,blocksize);
#endif
      SKIP;
      if(strcmp(label,blocklabel)!=0)
      {
        fseek(fd,blocksize,1);
        blocksize=0;
      }
    }
  }
  return(blocksize-8);
}

/* THE FOLLOWING METHODS WHERE TAKEN FROM THE readgadget.c FILE DISTRIBUTED BY
 KLAUS DOLAG */


/*---------------------- Basic routine to read data from a file ---------------*/
 template<typename PType>
 size_t GadgetFile<PType>::my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
 {
 size_t nread;
 
 if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
 {
 printf("I/O error (fread) !\n");
 exit(3);
 }
 return nread;
 }


template<typename PType>
bool GadgetFile<PType>::FileExists(string strFilename) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;
  
  // Attempt to get the file attributes
  intStat = stat(strFilename.c_str(),&stFileInfo);
  if(intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }
  
  return(blnReturn);
}

template<typename PType>
void GadgetFile<PType>::readBlock(const char *blockname) {
  
  float umass = 1.0; // ????
  
  int n;
  int np_in_file=np_file_end-np_file_start+1;
  cout << string(blockname) << endl;
  if (string(blockname)=="POS") {
    cout << "Reading block "<< string(blockname) << endl;
    pvector *pos=new pvector[ntot];
    cout << "...memory allocated." << endl;
    n = read_gadget_float3((float*)pos,"POS ",fd);
    cout << "...data are in memory." << endl;
    for (int i=0; i<=np_in_file-1; i++) {
      //pv[i+np_file_start]->setPos(pos[i].x,pos[i].y,pos[i].z);
      p_data[i+np_file_start][0] = pos[i].x;
      p_data[i+np_file_start][1] = pos[i].y;
      p_data[i+np_file_start][2] = pos[i].z;
    }
    delete [] pos;
    return;
  }  else if (string(blockname)=="VEL") {
    cout << "Reading block "<< string(blockname) << endl;
    pvector *vel=new pvector[ntot];
    cout << "...memory allocated..." << endl;
    n = read_gadget_float3((float*)vel,"VEL ",fd);
    cout << "...data is in memory." << endl;
    for (int i=0; i<=np_in_file-1; i++) {
      //pv[i+np_file_start]->setVel(vel[i].x,vel[i].y,vel[i].z);
      p_data[i+np_file_start].Vel[0] = vel[i].x;
      p_data[i+np_file_start].Vel[1] = vel[i].y;
      p_data[i+np_file_start].Vel[2] = vel[i].z;
    }
    delete [] vel;
    return;
  } else if (string(blockname)=="MASS"){
    cout << "Reading block "<< string(blockname) << endl;
    int nmass=0;
    for (int i=0; i<=5; i++) {
      if (masstab[i]==0.0) {nmass=nmass+npart[i];}}
    float *mass=new float[nmass];
    cout << "...memory allocated..." << endl;
    n = read_gadget_float((float*)mass,"MASS",fd);
    cout << "...data is in memory." << endl;
    int icounter=0;
    for (int i=0; i<=np_in_file-1; i++) {
      int ips=p_data[i+np_file_start].type;
      if (masstab[ips] > 0.0) {
        //pv[i+np_file_start]->setMass(masstab[ips]*umass);
        p_data[i+np_file_start].Mass = masstab[ips]*umass;
      } else {
        //pv[i+np_file_start]->setMass(mass[icounter]*umass);
        p_data[i+np_file_start].Mass = mass[icounter]*umass;
        icounter++;
      }
    }
    delete [] mass;
    return;
  }else if (string(blockname)=="RHO"){
    cout << "Reading block "<< string(blockname) << endl;
    float *rho=new float[npart[0]];
    cout << "...memory allocated..." << endl;
    n = read_gadget_float((float*)rho,"RHO ",fd);
    cout << "...data is in memory." << endl;
    for (int i=0; i<=npart[0]-1; i++) {
      //pv[i+np_file_start]->setRho(rho[i]);
      p_data[i+np_file_start].Rho = rho[i];
    }
    delete [] rho;
    return;
  } else if (string(blockname)=="RHOD"){
    cout << "Reading block "<< string(blockname) << endl;
    float *rho=new float[npart[1]];
    cout << "...memory allocated..." << endl;
    n = read_gadget_float((float*)rho,"RHOD",fd);
    cout << "...data is in memory." << endl;
    int icount=0;
    for (int i=npart[0]; i<=npart[0]+npart[1]-1; i++) {
      //pv[i+np_file_start]->setRho(rho[icount]*umass);
      p_data[i+np_file_start].Rho = rho[icount]*umass;
      icount++;
    }
    delete [] rho;
    return;
  }   else if (string(blockname)=="TEMP"){
    cout << "Reading block "<< string(blockname) << endl;
    float *temp=new float[npart[0]];
    cout << "...memory allocated..." << endl;
    n = read_gadget_float((float*)temp,"TEMP",fd);
    cout << "...data is in memory." << endl;
    int icount=0;
    for (int i=0; i<=npart[0]-1; i++) {
      //pv[i+np_file_start]->setTemp(temp[icount]);
      p_data[i+np_file_start].Temp = temp[icount];
      icount++;
    }
    delete [] temp;
    return;
  }else{
    std::cerr << "*** No " << blockname << " in file(s) " << filebasename << " ***" << std::endl;
  }
  
  //
}

template<>
void GadgetFile<ParticleType<float> >::readBlock(const char *blockname) {
  
  float umass = 1.0; // ????
  
  int n;
  int np_in_file=np_file_end-np_file_start+1;
  cout << string(blockname) << endl;
  if (string(blockname)=="POS") {
    cout << "Reading block "<< string(blockname) << endl;
    pvector *pos=new pvector[ntot];
    cout << "...memory allocated." << endl;
    n = read_gadget_float3((float*)pos,"POS ",fd);
    cout << "...data are in memory." << endl;
    for (int i=0; i<=np_in_file-1; i++) {
      //pv[i+np_file_start]->setPos(pos[i].x,pos[i].y,pos[i].z);
      p_data[i+np_file_start][0] = pos[i].x;
      p_data[i+np_file_start][1] = pos[i].y;
      p_data[i+np_file_start][2] = pos[i].z;
    }
    delete [] pos;
    return;
  } else if (string(blockname)=="MASS"){
    cout << "Reading block "<< string(blockname) << endl;
    int nmass=0;
    for (int i=0; i<=5; i++) {
      if (masstab[i]==0.0) {nmass=nmass+npart[i];}}
    float *mass=new float[nmass];
    cout << "...memory allocated..." << endl;
    n = read_gadget_float((float*)mass,"MASS",fd);
    cout << "...data is in memory." << endl;
    int icounter=0;
    for (int i=0; i<=np_in_file-1; i++) {
      int ips=p_data[i+np_file_start].type;
      if (masstab[ips] > 0.0) {
        //pv[i+np_file_start]->setMass(masstab[ips]*umass);
        p_data[i+np_file_start].Mass = masstab[ips]*umass;
      } else {
        //pv[i+np_file_start]->setMass(mass[icounter]*umass);
        p_data[i+np_file_start].Mass = mass[icounter]*umass;
        icounter++;
      }
    }
    delete [] mass;
    return;
  }else{
    std::cerr << "*** You cannot load " << blockname << " from file(s) " << filebasename << " into a ParticleType type. ***" << std::endl;
  }
}


/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to read a 1D float array ---------------------*/
/*-------- int *data:     Pointer where the data are stored to ----------------*/
/*-------- char *label:   Identifyer for the datafield to read ----------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns length of dataarray ----------------------------------------*/
/*-----------------------------------------------------------------------------*/
template<typename PType>
int GadgetFile<PType>::read_gadget_float(float *data,const char *label,FILE *fd)
{
  int blocksize;
  
  blocksize = find_block(fd,label);
  if(blocksize <= 0)
  {
    printf("Block <%s> not fond !\n",label);
    exit(5);
  }
  else
  {
#ifdef MY_DEBUG
    printf("Reading %d bytes of data from <%s>...\n",blocksize,label);
#endif
    SKIP;
    my_fread(data,blocksize, 1, fd);
    swap_Nbyte((char*)data,blocksize/sizeof(float),4);
    SKIP;
  }
  return(blocksize/sizeof(float));
}

/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to read a 3D float array ---------------------*/
/*-------- int *data:     Pointer where the data are stored to ----------------*/
/*-------- char *label:   Identifyer for the datafield to read ----------------*/
/*-------- FILE *fd:      File handle -----------------------------------------*/
/*-------- returns length of dataarray ----------------------------------------*/
/*-----------------------------------------------------------------------------*/
 template<typename PType>
 int GadgetFile<PType>::read_gadget_float3(float *data,const char *label,FILE *fd)
 {
 int blocksize;
 
 blocksize = find_block(fd,label);
 if(blocksize <= 0)
 {
 printf("Block <%s> not fond !\n",label);
 exit(5);
 }
 else
 {
 #ifdef MY_DEBUG
 printf("Reding %d bytes of data from <%s>...\n",blocksize,label);
 #endif
 SKIP;
 my_fread(data,blocksize, 1, fd);
 swap_Nbyte((char*)data,blocksize/sizeof(float),4);
 SKIP;
 }
 return(blocksize/sizeof(float)/3);
 }


#endif
