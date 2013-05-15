/*
 * model.h
 *
 *  Created on: Jan 23, 2012
 *      Author: mpetkova
 */

#ifndef MODEL_H_
#define MODEL_H_

#include <singlelens.h>
#include <multiplane.h>
#include <source.h>
#include "parameters.h"

/**
 * \brief A class that allocates the lens, the source, and the cosmological models.
 * There are several constructors, which allow for different types of initializations.
 *
 * The template class will then be created from the main() function as:
 * e.g. Model<MultiLens,SourceUniform>* model = new Model<MultiLens,SourceUniform>(paramfile,&seed);
 *
 * The default values of the Model template are AnaLens and SourceUnifrom, which can de allocated as
 * Model<>* model = new Model<>(paramfile);
 * equivalent to
 * Model<AnaLens,SourceUniform>* model = new Model<AnaLens,SourceUniform>(paramfile);
 *
 * the parameters correspond to the ones that the lens need. For example
 * MultiPlane == paramfile + seed
 * AnaLens == paramfile
 * etc.
 *
 */
template <class L=SingleLens, class S=SourceUniform> class Model{
public:

	L* lens;
	S* source;
	CosmoHndl cosmo;
	InputParams* params;

	///For the MultiPlane one
	Model(std::string paramfile,long* seed){
		params = new InputParams(paramfile);

		cosmo = new COSMOLOGY();
		lens = new L(*params,seed);
		source = new S(*params);

		setInternal();
	};
	/// For the AnaLens one
	Model(std::string paramfile){
		params = new InputParams(paramfile);

		cosmo = new COSMOLOGY();
		lens = new L(*params);
		source = new S(*params);

		setInternal();
	};
	///Others
	Model(long* seed){
		params = NULL;

		cosmo = new COSMOLOGY();
		lens = new L();
		source = new S();

		setInternal();
	};
	Model(){
		params = NULL;
		cosmo = new COSMOLOGY();
		lens = new L();
		source = new S();

		setInternal();
	};
	~Model(){
		if(params)
			delete params;
		delete lens;
		delete source;
		delete cosmo;
	};

	double getZsource(){return source->getZ();};
	double getZlens(){return lens->getZlens();};

    void RandomizeModel(double r_source_physical,long *seed
    		,bool tables,bool randomize_host_z=true,bool randomize_source_z=true,bool in_radians=false);
	
	inline void randomize(double step, long* seed)
	{
		//lens->randomize(step, seed);
		source->randomize(step, seed);
	}
	
	// write parameters
	void setParameters(Parameters& p)
	{
		//lens->setParameters(p);
		source->setParameters(p);
	}
	
	// read parameters
	void getParameters(Parameters& p)
	{
		//lens->getParameters(p);
		source->getParameters(p);
	}
	
private:

    void setInternal(){
    	lens->setInternalParams(cosmo,source);
    	//source->setDlDs(cosmo->angDist(0,lens->getZlens()) / cosmo->angDist(0,source->getZ()));
    };

    //void change_redshifts(TreeHndl i_tree,TreeHndl s_tree,double z_source,double z_lens);
};

/** \ingroup ChangeLens
* \brief routines for randomizing the lens.  How the lens is randomized is specified in the specific
* derived lens class that was used to construct the model.
 *
 */
template<class L,class S> void Model<L,S>::RandomizeModel(
		double r_source_phys
		,long *seed
		,bool tables
		,bool randomize_host_z
		,bool randomize_source_z
		,bool in_radians
		){
	double *zlTable,*zsTable;
	int n,i,NzTable;
	std::ifstream file;
	char *filename;

	if(tables){

		std::cout << "reading lens distribution tables" << std::endl;

		// read in Einstein radius projected onto the source plane

		filename = "GalaxyData/z_table.txt";
		file.open(filename);

		if(!file){
			std::cout << "Can't open file " << filename << std::endl;
			exit(1);
		}

		file >> NzTable;

		zsTable=new double[NzTable];
		zlTable=new double[NzTable];

		for(n=0;n<NzTable;++n) file >> zsTable[n] >> zlTable[n];

		file.close();

		// choose random set of redshifts
		double zlens, zsource;
		do{
			zlens = Utilities::RandomFromTable(zlTable,NzTable,seed);
			zsource = Utilities::RandomFromTable(zsTable,NzTable,seed);
		}while(zsource < zlens);

		delete[] zsTable;
		delete[] zlTable;

		if(randomize_source_z) source->setZ(zsource);
		if(randomize_host_z) lens->setZlens(cosmo,zlens,zsource);

		lens->RandomizeSigma(seed,tables);
	}

	// This randomizes the halos if they are not read from an external source
	setInternal();

	// This need to be done after source->DlDs has been set in setInternal()
	if(in_radians)
		source->setRadius(r_source_phys/cosmo->angDist(0,source->getZ()));
	else
		source->setRadius(r_source_phys*source->getDlDs());

	lens->RandomizeHost(seed,tables);

	return ;
}

#endif /* MODEL_H_ */
