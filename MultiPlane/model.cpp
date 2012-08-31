/*
 * model.cpp
 *
 *  Created on: Jan 23, 2012
 *      Author: mpetkova
 */

#include <model.h>
#include <source_models.h>
#include <sstream>
#include <string>
#include <fstream>

using namespace std;

/// TODO: MARGARITA/BEN the constructor has to be done in a smart way, so that it knows what kind of Lens and Source to create
/*
 * This creates a model and constructs a Lens, a Source, and a Cosmology.
 * The Type of the lens is determined by multi_lens, the type of the source is read in from the parameter file
 * and the the cosmology is the standard one
 */
Model::Model(std::string paramfile  /// parameterfile
		,long *my_seed		/// the seed, set in the main()
		,bool multi_lens	/// if true, then create a multi lens plane, if false -- stick to an analytic lens and single plane
		){

	readParamfile(paramfile);

	if(multi_lens == true){
		lens = new MultiLens(paramfile,my_seed);
		multilens = static_cast<MultiLens*>(lens);
	}
	else{
		lens = new AnaLens(paramfile);
		analens = static_cast<AnaLens*>(lens);
	}

	switch(sb_type){
	case Uniform:
		source = new SourceUniform(paramfile);
		break;
	case Gaussian:
		source = new SourceGaussian(paramfile);
		break;
	case BLR_Disk:
		source = new SourceBLRDisk(paramfile);
		break;
	case BLR_Sph1:
		source = new SourceBLRSph1(paramfile);
		break;
	case BLR_Sph2:
		source = new SourceBLRSph2(paramfile);
		break;
	default:
		std::cout << "!!!unrecongizable source type!!!" << std::endl;
		break;

	}

	cosmo = new COSMOLOGY();

	/// sets some source and lens paramaters that are interdependent
	/// e.g. source->DlDs depens on source and lens
	setInternal();
}

Model::~Model(){
	delete lens;
	delete source;
	delete cosmo;
}

void Model::readParamfile(std::string filename){
	  string label, rlabel, rvalue;
	  void *addr;
	  stringstream ss;
	  int i ,n, id = 1;
	  int mysbtype;
	  char dummy[100];
	  string escape = "#";
	  int flag;

	  addr = &sb_type;
	  label = "SourceSBType";

	  ifstream file_in(filename.c_str());
	  if(!file_in){
	    cout << "Can't open file " << filename << endl;
	    exit(1);
	  }

	  // output file
	  while(!file_in.eof()){
		  file_in >> rlabel >> rvalue;
		  file_in.getline(dummy,100);

		  if(rlabel[0] == escape[0])
			  continue;

		  flag = 0;

		  if(rlabel == label){

			  flag = 1;
			  ss << rvalue;

			  ss >> mysbtype;
			  *((int *)addr) = mysbtype;

			  id = -1;
		  }
	  }


	  if(id >= 0){
		  ERROR_MESSAGE();
		  cout << "parameter " << label << " needs to be set in the parameter file " << filename << endl;
		  exit(1);
	  }

	  file_in.close();

}

void Model::setInternal(){
	lens->setInternalParams(cosmo,source);
	source->setDlDs(cosmo->angDist(0,lens->getZlens()) / cosmo->angDist(0,source->getZ()));
}

