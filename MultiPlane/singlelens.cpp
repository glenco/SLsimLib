/*
 * singlelens.h
 *
 *  Created on: Aug 8, 2012
 *      Author: mpetkova
 */

#include <slsimlib.h>
#include <sstream>
#include <string>

using namespace std;

/*
 * \ingroup Constructor
 * allocates space for the halo trees and the inout lens, if there is any
 */
SingleLens::SingleLens(string filename,long *my_seed) : Lens(){
	readParamfile(filename);

	charge = 4*pi*Grav*mass_scale;

	seed = my_seed;

	/*if(internal_profile == NFW && tables_set != true){
		make_tables_nfw();
		tables_set = true;
	}
	if(internal_profile == PseudoNFW && tables_set != true){
		make_tables_pseudonfw(pnfw_beta);
		tables_set = true;
	}*/
}

SingleLens::~SingleLens(){
	delete halo_tree;
	delete halodata;

	/*if(tables_set == true){
		if(internal_profile == NFW) delete_tables_nfw();
		if(internal_profile == PseudoNFW) delete_tables_pseudonfw();
	}*/
}


void SingleLens::readParamfile(string filename){
      const int MAXPARAM = 50;
	  string label[MAXPARAM], rlabel, rvalue;
	  void *addr[MAXPARAM];
	  int id[MAXPARAM];
	  stringstream ss;
	  int i ,n;
	  int myint;
	  double mydouble;
	  string mystring;
	  char dummy[100];
	  string escape = "#";
	  int flag;

	  n = 0;

	  /// id[] = 2 = string, 1 = int, 0 = double
	  addr[n] = &outputfile;
	  id[n] = 2;
	  label[n++] = "outputfile";

	  addr[n] = &mass_scale;
	  id[n] = 0;
	  label[n++] = "mass_scale";

	  addr[n] = &internal_profile;
	  id[n] = 1;
	  label[n++] = "internal_profile";

	  addr[n] = &pw_beta;
	  id[n] = 0;
	  label[n++] = "internal_slope_pw";

	  addr[n] = &pnfw_beta;
	  id[n] = 0;
	  label[n++] = "internal_slope_pnfw";

	  addr[n] = &mass;
	  id[n] = 0;
	  label[n++] = "mass";

	  addr[n] = &zlens;
	  id[n] = 0;
	  label[n++] = "z_lens";

	  cout << "Single lens: reading from " << filename << endl;

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

		  for(i = 0; i < n; i++){
			  if(rlabel == label[i]){

				  flag = 1;
				  ss << rvalue;

				  switch(id[i]){
				  case 0:
					  ss >> mydouble;
					  *((double *)addr[i]) = mydouble;
					  break;
				  case 1:
					  ss >> myint;
					  *((int *)addr[i]) = myint;
					  break;
				  case 2:
					  ss >> mystring;
					  *((string *)addr[i]) = mystring;
					  break;
				  }

				  ss.clear();
				  ss.str(string());

				  id[i] = -1;
			  }
		  }
	  }

	  for(i = 0; i < n; i++){
		  if(id[i] >= 0 && addr[i] != &pw_beta && addr[i] != &pnfw_beta){
			  ERROR_MESSAGE();
			  cout << "parameter " << label[i] << " needs to be set!" << endl;
			  exit(0);
		  }

		  if(id[i] >= 0 && addr[i] == &pw_beta){
			  pw_beta = -1.0;
		  }
		  if(id[i] >= 0 && addr[i] == &pnfw_beta){
			  pnfw_beta = 2.0;
		  }
	  }

	  file_in.close();

	  if(pw_beta >= 0){
		  ERROR_MESSAGE();
		  cout << "Internal slope >=0 not possible." << endl;
		  exit(1);
	  }

	  if(pnfw_beta <= 0){
		  ERROR_MESSAGE();
		  cout << "Internal slope <=0 not possible." << endl;
		  exit(1);
	  }

	  if(pnfw_beta / floor(pnfw_beta) > 1.0){
		  ERROR_MESSAGE();
		  cout << "Internal slope needs to be a whole number." << endl;
		  exit(1);
	  }

	  printSingleLens();
}

void SingleLens::buildHaloTrees(
		CosmoHndl cosmo /// the cosmology
		,double zsource /// the source resdhift
		){

	halodata = new HaloData(cosmo,mass,zlens);

	switch(internal_profile){
	case PowerLaw:
		halo_tree = new ForceTreePowerLaw(pw_beta,&halodata->pos[0],halodata->Nhalos
				,halodata->halos,true,halodata->kappa_background);
		break;
	case NFW:
		halo_tree = new ForceTreeNFW(&halodata->pos[0],halodata->Nhalos
				,halodata->halos,true,halodata->kappa_background);
		break;
	case PseudoNFW:
		halo_tree = new ForceTreePseudoNFW(pnfw_beta,&halodata->pos[0],halodata->Nhalos
				,halodata->halos,true,halodata->kappa_background);
		break;
	default:
		cout << "There is no such case for the halo trees" << endl;
		ERROR_MESSAGE();
		exit(1);
		break;
	}

	/*
	switch(internal_profile){
	case PowerLaw:
		halo_tree = new QuadTreePowerLaw(pw_beta,&halodata->pos[0],halodata->Nhalos
								,halodata->halos,halodata->kappa_background);
		break;
	case NFW:
		halo_tree = new QuadTreeNFW(&halodata->pos[0],halodata->Nhalos
				,halodata->halos,halodata->kappa_background);
		break;
	case PseudoNFW:
		halo_tree = new QuadTreePseudoNFW(pnfw_beta,&halodata->pos[0],halodata->Nhalos
				,halodata->halos,halodata->kappa_background);
		break;
	default:
		cout << "There is no such case for the halo trees" << endl;
		ERROR_MESSAGE();
		exit(1);
		break;
	}
	*/
}

void SingleLens::rayshooterInternal(unsigned long Npoints, Point *i_points, bool kappa_off, double zsource){
	/* i_points need to be already linked to s_points */
	double kappa;
    double alpha[2], gamma[3];
    unsigned long i;

	for(i = 0; i< Npoints; i++){

		// do stars with tree code
		halo_tree->force2D(i_points[i].x,alpha,&kappa,gamma,kappa_off);

		i_points[i].image->x[0] += alpha[0];
		i_points[i].image->x[1] += alpha[1];

		if(!kappa_off){
			i_points[i].kappa = kappa;
			i_points[i].gamma[0] = gamma[0];
			i_points[i].gamma[1] = gamma[1];
			i_points[i].gamma[2] = 0.0;
		}

		// final operations to get the inverse magnification
    	i_points[i].invmag = (1-i_points[i].kappa)*(1-i_points[i].kappa)
  	    				- i_points[i].gamma[0]*i_points[i].gamma[0] - i_points[i].gamma[1]*i_points[i].gamma[1];
	}
}

void SingleLens::setInternalParams(CosmoHndl cosmo, SourceHndl source){
	buildHaloTrees(cosmo,source->getZ());
}


double SingleLens::getZlens(){
	return zlens;
}

void SingleLens::setZlens(double z){
	zlens = z;
}

void SingleLens::printSingleLens(){
	cout << endl << "outputfile " << outputfile << endl;

		cout << endl << "**single lens model**" << endl;

		cout << "mass scale " << mass_scale << endl;

		cout << "mass " << mass << endl;

		cout << "internal profile type " << internal_profile << endl;
		switch(internal_profile){
		case PowerLaw:
			cout << "  Power law internal profile " << endl;
			cout << "  slope: " << pw_beta << endl;
			break;
		case NFW:
			cout << "  NFW internal profile " << endl;
			break;
		case PseudoNFW:
			cout << "  Pseudo NFW internal profile " << endl;
			cout << "  slope: " << pnfw_beta << endl;
			break;
		}

		cout << endl;
}

void SingleLens::RandomizeHost(long *seed,bool tables){};
void SingleLens::RandomizeSigma(long *seed,bool tables){};
