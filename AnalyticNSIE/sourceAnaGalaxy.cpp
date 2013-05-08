/*
 * sourceAnaLens.cpp
 *
 *  Created on: Aug 13, 2012
 *      Author: bmetcalf
 */
#include "slsimlib.h"

// TODO: set `mag_limit` and `band` to default values in all constructors

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
	
	galaxies.push_back(OverzierSource(mag,BtoT,Reff,Rh,PA,inclination,0,my_z,my_theta));
}
/** Constructor for passing in a pointer to the galaxy model or a list of galaxies instead of constructing it internally.
*   Useful when there is a list of pre-allocated sources.  The redshifts and sky positions need to be set separately.
*/
MultiSourceAnaGalaxy::MultiSourceAnaGalaxy(
		OverzierSource *my_galaxy
		): Source(),index(0){
	
	galaxies.push_back(*my_galaxy);
}
/// Constructor for importing from data file.
MultiSourceAnaGalaxy::MultiSourceAnaGalaxy(
		InputParams& params   /// Input data file for galaxies
		): Source(),index(0){

	if(!params.get("input_galaxy_file",input_gal_file)){
		std::cout << "ERROR: input_galaxy_file not found in parameter file " << params.filename() << std::endl;
		exit(1);
	}

	assignParams(params);

	std::cout << "Constructing SourceAnaGalaxy" << std::endl;

	readDataFile();
	index = 0;
}

MultiSourceAnaGalaxy::~MultiSourceAnaGalaxy()
{
}

/// read in galaxies from a Millennium simulation file
void MultiSourceAnaGalaxy::readDataFile(){

	char c='0';
	//int type;
	//long galid,haloid;
	/*double cx,cy,cz,ra,dec,z_geo ,z_app ,dlum ,vlos ,incl
	,acs435 ,acs606 ,acs775,acs850,acs435_bulge,acs606_bulge,acs775_bulge,acs850_bulge
	,mvir,rvir,vmax,stellarmass,bulgemass,stellardiskradius,bulgesize
	,inclination,pa,angdist,diskradius_arcsec,bulgesize_arcsec;*/
	unsigned long i,j;

	unsigned long GalID,HaloID;
	double ra,dec,z_cosm,z_app,Dlum,inclination,pa,Rh,Ref,SDSS_u,SDSS_g,SDSS_r,SDSS_i,SDSS_z
	,J_band,H_band,Ks_band,i1,i2,SDSS_u_Bulge,SDSS_g_Bulge,SDSS_r_Bulge,SDSS_i_Bulge,SDSS_z_Bulge
	,J_band_Bulge,H_band_Bulge,Ks_band_Bulge,i1_Bulge,i2_Bulge;

	std::ifstream file_in(input_gal_file.c_str());
	if(!file_in){
		std::cout << "Can't open file " << input_gal_file << std::endl;
		exit(1);
	}

	std::cout << "Reading from galaxy data file " << input_gal_file.c_str() << std::endl;
	//file_in >> Ngalaxies;
	//std::cout << "Number of source galaxies: " << Ngalaxies << std::endl;

	//TODO Using a vector of pointer to OverzierSource is inefficient for mem.  check that a default copy would work for adding galaxies
	i=0;
	while(file_in.peek() == '#'){
		file_in.ignore(10000,'\n');
		++i;
	}
	std::cout << "skipped "<< i << " comment lines in " << input_gal_file << std::endl;

	double theta[2] = {0.0,0.0};

	const int ncolumns = 31;

	void *addr[ncolumns];
	addr[0] = &GalID;
	addr[1] = &HaloID;
	addr[2] = &ra;
	addr[3] = &dec;
	addr[4] = &z_cosm;
	addr[5] = &z_app;
	addr[6] = &Dlum;
	addr[7] = &inclination;
	addr[8] = &pa;
	addr[9] = &Rh;
	addr[10] = &Ref;
	addr[11] = &SDSS_u;
	addr[12] = &SDSS_g;
	addr[13] = &SDSS_r;
	addr[14] = &SDSS_i;
	addr[15] = &SDSS_z;
	addr[16] = &J_band;
	addr[17] = &H_band;
	addr[18] = &Ks_band;
	addr[19] = &i1;
	addr[20] = &i2;
	addr[21] = &SDSS_u_Bulge;
	addr[22] = &SDSS_g_Bulge;
	addr[23] = &SDSS_r_Bulge;
	addr[24] = &SDSS_i_Bulge;
	addr[25] = &SDSS_z_Bulge;
	addr[26] = &J_band_Bulge;
	addr[27] = &H_band_Bulge;
	addr[28] = &Ks_band_Bulge;
	addr[29] = &i1_Bulge;
	addr[30] = &i2_Bulge;

	unsigned long myint;
	double mydouble;
	std::string myline;
	std::string strg;
	std::string f=",";
	std::stringstream buffer;
	size_t length;

	double mag,mag_bulge;

	// read in data
	for(i=0,j=0 ; ; ++i){
		myline.clear();
		getline(file_in,myline);

		if(myline[0] == '#')
			break;

		for(int l=0;l<ncolumns; l++){
			int pos = myline.find(f);
			strg.assign(myline,0,pos);
			buffer << strg;
			if(l == 0 || l == 1 ){
				buffer >> myint;
				*((unsigned long *)addr[l]) = myint;
			}
			else{
				buffer >> mydouble;
				*((double *)addr[l]) = mydouble;
			}
			myline.erase(0,pos+1);
			strg.clear();
			buffer.clear();
			buffer.str(std::string());
		}

		// read a line of data
		/*file_in >> galid >> c >> haloid >> c >> cx >> c >> cy >> c >> cz >> c >> ra >> c >> dec >> c >> z_geo >> c >> z_app
		>> c >> dlum >> c >> vlos >> c >> incl
		>> c >> acs435 >> c >> acs606 >> c >> acs775 >> c >> acs850
		>> c >> acs435_bulge >> c >> acs606_bulge >> c >> acs775_bulge >> c >> acs850_bulge
		>> c >> type >> c >> mvir >> c >> rvir >> c >> vmax >> c >> stellarmass >> c >> bulgemass
		>> c >> stellardiskradius >> c >> bulgesize
		>> c >> inclination >> c >> pa >> c >> angdist >> c >> diskradius_arcsec >> c >> bulgesize_arcsec >> c;
		*/
		/*file_in >> GalID  >> HaloID  >> ra  >> dec  >> z_cosm  >> z_app  >> Dlum  >> inclination
				 >> pa  >> Rh  >> Ref  >> SDSS_u  >> SDSS_g  >> SDSS_r  >> SDSS_i  >> SDSS_z
				 >> J_band  >> H_band  >> Ks_band  >> i1  >> i2  >> SDSS_u_Bulge  >> SDSS_g_Bulge  >> SDSS_r_Bulge
				 >> SDSS_i_Bulge  >> SDSS_z_Bulge  >> J_band_Bulge  >> H_band_Bulge  >> Ks_band_Bulge  >> i1_Bulge
				 >> i2_Bulge ;

		file_in >> GalID >> c >> HaloID >> c >> ra >> c >> dec >> c >> z_cosm >> c >> z_app >> c >> Dlum >> c >> inclination
		>> c >> pa >> c >> Rh >> c >> Ref >> c >> SDSS_u >> c >> SDSS_g >> c >> SDSS_r >> c >> SDSS_i >> c >> SDSS_z
		>> c >> J_band >> c >> H_band >> c >> Ks_band >> c >> i1 >> c >> i2 >> c >> SDSS_u_Bulge >> c >> SDSS_g_Bulge >> c >> SDSS_r_Bulge
		>> c >> SDSS_i_Bulge >> c >> SDSS_z_Bulge >> c >> J_band_Bulge >> c >> H_band_Bulge >> c >> Ks_band_Bulge >> c >> i1_Bulge
		>> c >> i2_Bulge >> c;
*/

		addr[11] = &SDSS_u;
		addr[12] = &SDSS_g;
		addr[13] = &SDSS_r;
		addr[14] = &SDSS_i;
		addr[15] = &SDSS_z;
		addr[16] = &J_band;
		addr[17] = &H_band;
		addr[18] = &Ks_band;
		addr[19] = &i1;
		addr[20] = &i2;

		switch (band){
		case SDSS_U:
			mag = SDSS_u;
			mag_bulge = SDSS_u_Bulge;
			break;
		case SDSS_G:
			mag = SDSS_g;
			mag_bulge = SDSS_g_Bulge;
			break;
		case SDSS_R:
			mag = SDSS_r;
			mag_bulge = SDSS_r_Bulge;
			break;
		case SDSS_I:
			mag = SDSS_i;
			mag_bulge = SDSS_i_Bulge;
			break;
		case SDSS_Z:
			mag = SDSS_z;
			mag_bulge = SDSS_z_Bulge;
			break;
		case J:
			mag = J_band;
			mag_bulge = J_band_Bulge;
			break;
		case H:
			mag = H_band;
			mag_bulge = H_band_Bulge;
			break;
		case Ks:
			mag = Ks_band;
			mag_bulge = Ks_band_Bulge;
			break;
		/*case i1:
			mag = i1;
			mag_bulge = i1_Bulge;
			break;
		case i2:
			mag = i2;
			mag_bulge = i2_Bulge;
			break;*/
		default:
			std::cout << "Requested band is not an available option." << std::endl;
			ERROR_MESSAGE();
			exit(1);
		}
		if(mag < mag_limit){
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
      
      if(j == 0){
        rangex[0] = rangex[1] = theta[0];
        rangey[0] = rangey[1] = theta[1];
      }else{
        if(theta[0] < rangex[0]) rangex[0] = theta[0];
        if(theta[0] > rangex[1]) rangex[1] = theta[0];
        if(theta[1] < rangey[0]) rangey[0] = theta[1];
        if(theta[1] > rangey[1]) rangey[1] = theta[1];
      }
      
			/***************************/
			galaxies.push_back(
					OverzierSource(mag,pow(10,-(mag_bulge-mag)/2.5),Ref,Rh
							,pa,inclination,HaloID,z_cosm,theta)
			);

			galaxies.back().setUMag(SDSS_u);
			galaxies.back().setGMag(SDSS_g);
			galaxies.back().setRMag(SDSS_r);
			galaxies.back().setIMag(SDSS_i);
			galaxies.back().setZMag(SDSS_z);
			galaxies.back().setJMag(J_band);
			galaxies.back().setHMag(H_band);
			galaxies.back().setKMag(Ks_band);

			//std::cout << "z:" << z_cosm << " mag " << SDSS_u << " Bulge to total " << pow(10,-(SDSS_u_Bulge-SDSS_u)/2.5)
			//		<< " bulge size arcsec " << Ref  << " disk size arcsec " << pa << " position angle " << pa << " inclination " << inclination
			//		<< " theta = " << theta[0] << " " << theta[1] << std::endl;

			++j;
		}
	}

	file_in.close();

	std::cout << galaxies.size() << " galaxies read in."<< std::endl;
	std::cout << "angular range x = " << rangex[0] << " to " << rangex[1] << "   y = " << rangey[0] << " to " << rangey[1]
  << std::endl;

	return;
}

void MultiSourceAnaGalaxy::assignParams(InputParams& params){
	if(!params.get("input_galaxy_file",input_gal_file)){
		  ERROR_MESSAGE();
		  std::cout << "parameter input_galaxy_file needs to be set in the parameter file "
				  << params.filename() << std::endl;
		  exit(0);
	  }
	if(!params.get("source_band",band)){
		std::cout << "ERROR: Must assign source_band in parameter file " << params.filename() << std::endl;
		std::cout << "Could be that specified band is not available " << std::endl;
		exit(1);
	}
	if(!params.get("source_mag_limit",mag_limit)){
		std::cout << "ERROR: Must assign source_mag_limit in parameter file " << params.filename() << std::endl;
		exit(1);
	}

	if(!params.get("source_sb_limit",sb_limit))
		setSBlimit_magarcsec(30.);
	else
		sb_limit = pow(10,-0.4*(48.6+sb_limit))*pow(180*60*60/pi,2)/hplanck;
		/*{
		std::cout << "ERROR: Must assign source_sb_limit in parameter file " << params.filename() << std::endl;
		exit(1);*/
}
/** \brief Artificially increase the number of sources to increase the incidence of strong lenses.
 *
 * Sources above the limiting redshift and below the limiting magnitude are copied and given random
 * position angles, inclinations and positions within the field of view.  The field of view is determined
 * directly from the range of already existing source positions.  The field is assumed to be rectangular.
 */
void MultiSourceAnaGalaxy::multiplier(
		double z                /// limiting redshift, only sources above this redshift are copied
		,double mag_cut         /// limiting magnitude, only sources with magnitudes below this limit will be copied
		,int multiplicity       /// the number of times each of these sources should be multiplied
		,long *seed    /// random number seed
		){
	unsigned long Nold = galaxies.size(),NtoAdd=0;
	double x1[2],x2[2],theta[2];

	for(unsigned long i=0;i<Nold;++i){
		if(galaxies[i].getX()[0] < x1[0]) x1[0] = galaxies[i].getX()[0];
		if(galaxies[i].getX()[0] > x2[0]) x2[0] = galaxies[i].getX()[0];
		if(galaxies[i].getX()[1] < x1[1]) x1[1] = galaxies[i].getX()[1];
		if(galaxies[i].getX()[1] > x2[1]) x2[1] = galaxies[i].getX()[1];

		if(galaxies[i].getZ() > z && galaxies[i].getMag() < mag_cut) ++NtoAdd;
	}

	NtoAdd = 0;
	for(unsigned long i=0;i<Nold;++i){
		if(galaxies[i].getZ() > z && galaxies[i].getMag() < mag_cut){

			for(int j=0;j<multiplicity;++j){
				theta[0] = x1[0] + (x2[0] - x1[0])*ran2(seed);
				theta[1] = x1[1] + (x2[1] - x1[1])*ran2(seed);

				galaxies.push_back(OverzierSource(galaxies[i].getMag(),galaxies[i].getBtoT(),galaxies[i].getReff()
					,galaxies[i].getRh(),ran2(seed)*pi,ran2(seed)*2*pi
					,Nold+NtoAdd,galaxies[i].getZ(),theta));

				++NtoAdd;
			}
		}
	}

}

/// Sort the sources by redshift in assending order
void MultiSourceAnaGalaxy::sortInRedshift(){
    std::sort(galaxies.begin(),galaxies.end(),redshiftcompare);
}
// used in MultiSourceAnaGalaxy::sortInRedshift()
bool redshiftcompare(OverzierSource s1,OverzierSource s2){
    return (s1.getZ() < s2.getZ());
}
/// Sort the sources by magnitude in assending order
void MultiSourceAnaGalaxy::sortInMag(){
    std::sort(galaxies.begin(),galaxies.end(),magcompare);
}
// used in MultiSourceAnaGalaxy::sortInRedshift()
bool magcompare(OverzierSource s1,OverzierSource s2){
    return (s1.getMag() < s2.getMag());
}
/// Print info on current source parameters
void MultiSourceAnaGalaxy::printSource(){
	std::cout << "Overzier Galaxy Model" << std::endl;
	galaxies[index].printSource();
}
