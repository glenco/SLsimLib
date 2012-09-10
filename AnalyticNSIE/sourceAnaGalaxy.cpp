/*
 * sourceAnaLens.cpp
 *
 *  Created on: Aug 13, 2012
 *      Author: bmetcalf
 */
#include <slsimlib.h>
#include <sstream>

/// Source model for a single analytic galaxy model.
MultiSourceAnaGalaxy::MultiSourceAnaGalaxy(
		double mag              /// Total magnitude
		,double BtoT            /// Bulge to total ratio
		,double Reff            /// Bulge half light radius (arcs)
		,double Rh              /// disk scale hight (arcs)
		,double PA              /// Position angle (radians)
		,double inclination     /// inclination of disk (radians)
		,double my_z               /// redshift of source
		,double *my_theta          /// position on the sky
		): Source(),index(0){

	source_sb_type = MultiAnaSource;
	mem_allocated = true;
	galaxies.push_back(new OverGalaxy(mag,BtoT,Reff,Rh,PA,inclination,my_z,my_theta));
}
/** Constructor for passing in a pointer to the galaxy model or a list of galaxies instead of constructing it internally.
*   Useful when there is a list of pre-allocated sources.  The redshifts and sky positions need to be set separately.
*/
MultiSourceAnaGalaxy::MultiSourceAnaGalaxy(
		OverGalaxy *my_galaxy
		): Source(),index(0){

	source_sb_type = MultiAnaSource;
	mem_allocated = false;
	galaxies.push_back(my_galaxy);
}
/// Constructor for importing from data file.
MultiSourceAnaGalaxy::MultiSourceAnaGalaxy(
		std::string filename   /// Input data file for galaxies
		,double my_mag_limit   /// Apparent Magnitude limit
		){

	source_sb_type = MultiAnaSource;
	readParamfile(filename);

	mem_allocated = true;
	std::cout << "Constructing SourceAnaGalaxy" << std::endl;

	readDataFile(input_gal_file,my_mag_limit);
}

MultiSourceAnaGalaxy::~MultiSourceAnaGalaxy(){
	if(mem_allocated)
		for(unsigned long i=0;i<galaxies.size();++i) delete galaxies[i];
}

/// read in galaxies from a Millennium simulation file
void MultiSourceAnaGalaxy::readDataFile(std::string input_gal_file,double mag_limit){

	char c='0';
	//int type;
	//long galid,haloid;
	/*double cx,cy,cz,ra,dec,z_geo ,z_app ,dlum ,vlos ,incl
	,acs435 ,acs606 ,acs775,acs850,acs435_bulge,acs606_bulge,acs775_bulge,acs850_bulge
	,mvir,rvir,vmax,stellarmass,bulgemass,stellardiskradius,bulgesize
	,inclination,pa,angdist,diskradius_arcsec,bulgesize_arcsec;*/
	unsigned long i,j,Ngalaxies;

	unsigned long GalID,HaloID;
	double ra,dec,z_cosm,z_app,Dlum,inclination,pa,Rh,Ref,SDSS_u,SDSS_g,SDSS_r,SDSS_i,SDSS_z
	,J,H,Ks,i1,i2,SDSS_u_Bulge,SDSS_g_Bulge,SDSS_r_Bulge,SDSS_i_Bulge,SDSS_z_Bulge
	,J_Bulge,H_Bulge,Ks_Bulge,i1_Bulge,i2_Bulge;

	std::ifstream file_in(input_gal_file.c_str());
	if(!file_in){
		std::cout << "Can't open file " << input_gal_file << std::endl;
		exit(1);
	}

	std::cout << "Reading from galaxy data file " << input_gal_file.c_str() << std::endl;
	//file_in >> Ngalaxies;
	//std::cout << "Number of source galaxies: " << Ngalaxies << std::endl;

	//TODO Using a vector of pointer to OverGalaxy is inefficient for mem.  check that a default copy would work for adding galaxies
	//galaxies.resize(Ngalaxies);// = new OverGalaxy[Ngalaxies];
	i=0;
	while(file_in.peek() == '#'){
		file_in.ignore(10000,'\n');
		++i;
	}
	std::cout << "skipped "<< i << " comment lines in " << input_gal_file << std::endl;

	mag_limit = 25;  //TODO find a better way of passing this in
	double theta[2] = {0.0,0.0};

	// read in data
	for(i=0,j=0 ; c != '#' ; ++i){

		// read a line of data
		/*file_in >> galid >> c >> haloid >> c >> cx >> c >> cy >> c >> cz >> c >> ra >> c >> dec >> c >> z_geo >> c >> z_app
		>> c >> dlum >> c >> vlos >> c >> incl
		>> c >> acs435 >> c >> acs606 >> c >> acs775 >> c >> acs850
		>> c >> acs435_bulge >> c >> acs606_bulge >> c >> acs775_bulge >> c >> acs850_bulge
		>> c >> type >> c >> mvir >> c >> rvir >> c >> vmax >> c >> stellarmass >> c >> bulgemass
		>> c >> stellardiskradius >> c >> bulgesize
		>> c >> inclination >> c >> pa >> c >> angdist >> c >> diskradius_arcsec >> c >> bulgesize_arcsec >> c;
		*/

		file_in >> GalID >> c >> HaloID >> c >> ra >> c >> dec >> c >> z_cosm >> c >> z_app >> c >> Dlum >> c >> inclination
		>> c >> pa >> c >> Rh >> c >> Ref >> c >> SDSS_u >> c >> SDSS_g >> c >> SDSS_r >> c >> SDSS_i >> c >> SDSS_z
		>> c >> J >> c >> H >> c >> Ks >> c >> i1 >> c >> i2 >> c >> SDSS_u_Bulge >> c >> SDSS_g_Bulge >> c >> SDSS_r_Bulge
		>> c >> SDSS_i_Bulge >> c >> SDSS_z_Bulge >> c >> J_Bulge >> c >> H_Bulge >> c >> Ks_Bulge >> c >> i1_Bulge
		>> c >> i2_Bulge >> c;  //TODO the GalID will miss the first digit using this method.  No other method stops at the end of file.

			//TODO  BEN this needs to be selected from the parameter file
		if(SDSS_r < mag_limit){
			/*
			std::cout << galid << c << haloid << c << cx << c << cy << c << cz << c << ra << c << dec << c << z_geo << c << z_app
			<< c << dlum << c << vlos << c << incl
			<< c << acs435 << c << acs606 << c << acs775 << c << acs850
			<< c << acs435_bulge << c << acs606_bulge << c << acs775_bulge << c << acs850_bulge
			<< c << type << c << mvir << c << rvir << c << vmax << c << stellarmass << c << bulgemass
			<< c << stellardiskradius << c << bulgesize
			<< c << inclination << c << pa << c << angdist << c << diskradius_arcsec << c << bulgesize_arcsec << std::endl;
			 */

			theta[0] = ra*pi/180;
			theta[1] = dec*pi/180;

			/***************************/
			galaxies.push_back(
					new OverGalaxy(SDSS_r,pow(10,-(SDSS_r_Bulge-SDSS_r)/2.5),Ref,Rh
							,pa,inclination,z_cosm,theta)
			);
			//std::cout << "z:" << z[i] << " mag " << acs435 << " Bulge to total " << pow(10,-(acs435_bulge-acs435)/2.5) << " bulge size arcsec " << bulgesize_arcsec
			//		<< " disk size arcsec " << diskradius_arcsec << " position angle " << pa << " inclination " << inclination
			//		<< " theta = " << theta[i][0] << " " << theta[i][1] << std::endl;

			++j;
		}
	}

	file_in.close();

	std::cout << galaxies.size() << " galaxies read in."<< std::endl;

	Ngalaxies = galaxies.size();

	return;
}
void MultiSourceAnaGalaxy::readParamfile(std::string filename){
      const int MAXPARAM = 50;
	  std::string label[MAXPARAM], rlabel, rvalue;
	  void *addr[MAXPARAM];
	  int id[MAXPARAM];
	  std::stringstream ss;
	  int i ,n;
	  int myint;
	  double mydouble;
	  std::string mystring;
	  char dummy[100];
	  std::string escape = "#";
	  int flag;

	  n = 0;

	  /// id[] = 2 = string, 1 = int, 0 = double

	  addr[n] = &input_gal_file;
	  id[n] = 2;
	  label[n++] = "input_galaxy_file";

	  std::cout << "MultiSourceAnaGalaxy reading from " << filename << std::endl;

	  std::ifstream file_in(filename.c_str());
	  if(!file_in){
	    std::cout << "Can't open file " << filename << std::endl;
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
					  *((std::string *)addr[i]) = mystring;
					  break;
				  }

				  ss.clear();
				  ss.str(std::string());

				  id[i] = -1;
			  }
		  }
	  }

	  for(i = 0; i < n; i++){
		  if(id[i] >= 0){
			  ERROR_MESSAGE();
			  std::cout << "parameter " << label[i] << " needs to be set in the parameter file "
					  << filename << std::endl;
			  exit(0);
		  }
	  }

	  file_in.close();

}

/// Print info on current source parameters
void MultiSourceAnaGalaxy::printSource(){
	std::cout << "Overzier Galaxy Model" << std::endl;
	galaxies[index]->print();
}
