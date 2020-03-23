/*
 * sourceAnaLens.cpp
 *
 *  Created on: Aug 13, 2012
 *      Author: bmetcalf
 */
#include "slsimlib.h"
#include "simpleTreeVec.h"

// TODO: set `mag_limit` and `band` to default values in all constructors

/// Source model for a single analytic galaxy model.
SourceMultiAnaGalaxy::SourceMultiAnaGalaxy(
		PosType mag              /// Total magnitude
		,PosType mag_bulge        /// magnitude of Bulge
		,PosType Reff            /// Bulge half light radius (arcs)
		,PosType Rdisk              /// disk scale hight (arcs)
		,PosType PA              /// Position angle (radians)
		,PosType inclination     /// inclination of disk (radians)
		,PosType my_z               /// redshift of source
		,PosType *my_theta          /// position on the sky
    ,Utilities::RandomNumbers_NR &ran
		): Source(),index(0){
	
	galaxies.push_back(SourceOverzierPlus(mag,mag_bulge,Reff,Rdisk,PA,inclination,0,my_z,my_theta,ran));
}
/** Constructor for passing in a pointer to the galaxy model or a list of galaxies instead of constructing it internally.
*   Useful when there is a list of pre-allocated sources.  The redshifts and sky positions need to be set separately.
*/
SourceMultiAnaGalaxy::SourceMultiAnaGalaxy(
		SourceOverzierPlus *my_galaxy
		): Source(),index(0){
	
	galaxies.push_back(*my_galaxy);
}
// Constructor for importing from data file.
//SourceMultiAnaGalaxy::SourceMultiAnaGalaxy(
//		InputParams& params   /// Input data file for galaxies
//    ,Utilities::RandomNumbers_NR &ran
//		): Source(),index(0){
//
//	if(!params.get("source_input_galaxy_file",input_gal_file)){
//		std::cout << "ERROR: source_input_galaxy_file not found in parameter file " << params.filename() << std::endl;
//		exit(1);
//	}
//
//	assignParams(params);
//
//	std::cout << "Constructing SourceAnaGalaxy" << std::endl;
//
//	readDataFileMillenn(ran);
//	index = 0;
//  searchtree = new TreeSimpleVec<SourceOverzierPlus>(galaxies.data(),galaxies.size(),1,2,true,SourceOverzierPlus::getx);
//  //searchtree = new TreeSimpleVec<SourceOverzierPlus>(galaxies.data(),galaxies.size(),1,3,true);
//}

SourceMultiAnaGalaxy::~SourceMultiAnaGalaxy()
{
  delete searchtree;
}

/// read in galaxies from a Millennium simulation file
void SourceMultiAnaGalaxy::readDataFileMillenn(Utilities::RandomNumbers_NR &ran){

	//int type;
	//long galid,haloid;
	/*PosType cx,cy,cz,ra,dec,z_geo ,z_app ,dlum ,vlos ,incl
	,acs435 ,acs606 ,acs775,acs850,acs435_bulge,acs606_bulge,acs775_bulge,acs850_bulge
	,mvir,rvir,vmax,stellarmass,bulgemass,stellardiskradius,bulgesize
	,inclination,pa,angdist,diskradius_arcsec,bulgesize_arcsec;*/
	unsigned long i,j;

	unsigned long GalID,HaloID;
	PosType ra,dec,z_cosm,z_app,Dlum,inclination,pa,Rdisk,Ref,SDSS_u,SDSS_g
  ,SDSS_r,SDSS_i,SDSS_z,J_band,H_band,Ks_band,i1,i2
  ,SDSS_u_Bulge,SDSS_g_Bulge,SDSS_r_Bulge,SDSS_i_Bulge,SDSS_z_Bulge
	,J_band_Bulge,H_band_Bulge,Ks_band_Bulge,i1_Bulge,i2_Bulge;

  std::ofstream color_cat;
  color_cat.open(input_gal_file + "color_catalog");
  
  color_cat << "# column 1 ID number" << std::endl;
  color_cat << "# column 2 redshift" << std::endl;
  color_cat << "# column 3 SDSS_u" << std::endl;
  color_cat << "# column 4 SDSS_g" << std::endl;
  color_cat << "# column 5 SDSS_r" << std::endl;
  color_cat << "# column 6 SDSS_i" << std::endl;
  color_cat << "# column 7 SDSS_z" << std::endl;
  color_cat << "# column 8 J_band" << std::endl;
  color_cat << "# column 9 H_band" << std::endl;
  color_cat << "# column 10 Ks_band" << std::endl;
  color_cat << "# column 11 i1" << std::endl;
  color_cat << "# column 12 i2" << std::endl;
  color_cat << "# column 13 SDSS_u_Bulge" << std::endl;
  color_cat << "# column 14 SDSS_g_Bulge" << std::endl;
  color_cat << "# column 15 SDSS_r_Bulge" << std::endl;
  color_cat << "# column 16 SDSS_i_Bulge" << std::endl;
  color_cat << "# column 17 SDSS_z_Bulge" << std::endl;
  color_cat << "# column 18 J_band_Bulge" << std::endl;
  color_cat << "# column 19 H_band_Bulge" << std::endl;
  color_cat << "# column 20 Ks_band_Bulge" << std::endl;
  color_cat << "# column 21 i1_Bulge" << std::endl;
  color_cat << "# column 22 i2_Bulge" << std::endl;

	std::ifstream file_in(input_gal_file.c_str());
	if(!file_in){
		std::cout << "Can't open file " << input_gal_file << std::endl;
    ERROR_MESSAGE();
    throw std::runtime_error(" Cannot open file.");
		exit(1);
	}
  
  std::string myline;
  size_t count = 0;
  while(getline(file_in,myline)){
    if(myline[0] != '#') ++count;
  }
  file_in.clear();
  file_in.seekg(0);
  galaxies.reserve(count);

	std::cout << "Reading from galaxy data file " << input_gal_file.c_str() << std::endl;
	//file_in >> Ngalaxies;
	//std::cout << "Number of source galaxies: " << Ngalaxies << std::endl;

	//TODO: Using a vector of pointer to OverzierSource is inefficient for mem.  check that a default copy would work for adding galaxies
	i=0;
	while(file_in.peek() == '#'){
		file_in.ignore(10000,'\n');
		++i;
	}
	std::cout << "   skipped "<< i << " comment lines in " << input_gal_file << std::endl;
  getline(file_in,myline); // skip title line


	PosType theta[2] = {0.0,0.0};

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
	addr[9] = &Rdisk;
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
	PosType myPosType;
	std::string strg;
	std::string f=",";
	std::stringstream buffer;

	//PosType mag,mag_bulge;

	// read in data
	for(i=0,j=0 ; ; ++i){
		myline.clear();
		getline(file_in,myline);

		if(myline[0] == '#')
			break;

    //std::cout << "myline " << myline << std::endl ;

    //std::cout << "read columns" << std::endl;
		for(int l=0;l<ncolumns; l++){
			long pos = myline.find(f);
			strg.assign(myline,0,pos);
			buffer << strg;
			if(l == 0 || l == 1 ){
        //std::cout << "myint " << myint ;
				buffer >> myint;
				*((unsigned long *)addr[l]) = myint;
			}
			else{
        //std::cout << "myPosType " << myPosType ;
				buffer >> myPosType;
				*((PosType *)addr[l]) = myPosType;
			}
			myline.erase(0,pos+1);
			strg.clear();
			buffer.clear();
			buffer.str(std::string());
		}

  
    color_cat << GalID << "  " << z_app << "  " << SDSS_u << "  " <<SDSS_g << "  " <<SDSS_r << "  " <<
    SDSS_i << "  " << SDSS_z << "  " << J_band << "  " << H_band << "  " << Ks_band <<
    "  " << i1 << "  " << i2 << "  " << SDSS_u_Bulge << "  " << SDSS_g_Bulge << "  " <<
    SDSS_r_Bulge << "  " << SDSS_i_Bulge << "  " << SDSS_z_Bulge << "  " << J_band_Bulge
    << "  " << H_band_Bulge << "  " << Ks_band_Bulge << "  " << i1_Bulge << "  " << i2_Bulge
    << std::endl;

      // converting from Millennium conventions
			theta[0] = -ra*PI/180;
			theta[1] = dec*PI/180;
      pa = (90 - pa)*PI/180;
      inclination *= PI/180;
      if(cos(inclination)< 0.5) inclination = acos(0.5);
    
      //std::cout << "did inclination" << std::endl;
    if(j == 0){
      rangex[0] = rangex[1] = theta[0];
      rangey[0] = rangey[1] = theta[1];
    }else{
      //std::cout << "rang" << std::endl;

      if(theta[0] < rangex[0]) rangex[0] = theta[0];
      if(theta[0] > rangex[1]) rangex[1] = theta[0];
      if(theta[1] < rangey[0]) rangey[0] = theta[1];
      if(theta[1] > rangey[1]) rangey[1] = theta[1];
    }
      
      //std::cout << "adding to galaxies" << std::endl;
			/***************************/
    
    SourceOverzierPlus galaxy(SDSS_i,SDSS_i_Bulge,Ref,Rdisk
                                ,pa,inclination,HaloID,z_cosm,theta,ran);
     //std::cout << "filling last galaxy" << std::endl;

    galaxy.setMag(SDSS_U,SDSS_u);
    galaxy.setMagBulge(SDSS_U,SDSS_u_Bulge);
    galaxy.setMag(SDSS_G,SDSS_g);
    galaxy.setMagBulge(SDSS_G,SDSS_g_Bulge);
    galaxy.setMag(SDSS_R,SDSS_r);
    galaxy.setMagBulge(SDSS_R,SDSS_r_Bulge);
    galaxy.setMag(SDSS_I,SDSS_i);
    galaxy.setMagBulge(SDSS_I,SDSS_i_Bulge);
    galaxy.setMag(SDSS_Z,SDSS_z);
    galaxy.setMagBulge(SDSS_Z,SDSS_z_Bulge);
    galaxy.setMag(J,J_band);
    galaxy.setMagBulge(J,J_band_Bulge);
    galaxy.setMag(H,H_band);
    galaxy.setMagBulge(H,H_band_Bulge);
    galaxy.setMag(Ks,Ks_band);
    galaxy.setMagBulge(Ks,Ks_band_Bulge);

    galaxy.setMag(Ks,Ks_band);
    galaxy.setMagBulge(Ks,Ks_band_Bulge);

      // The Euclid bands are not actually read in
    galaxy.setMag(EUC_VIS,SDSS_i);
    galaxy.setMagBulge(EUC_VIS,SDSS_i_Bulge);
    galaxy.setMag(EUC_H,H_band);
    galaxy.setMagBulge(EUC_H,H_band_Bulge);
    galaxy.setMag(EUC_J,J_band);
    galaxy.setMagBulge(EUC_J,J_band_Bulge);

    galaxy.changeBand(band);
      
    if(galaxy.getMag() < mag_limit ){
      galaxies.push_back(galaxy);
      assert(galaxies.back().getMag(band) == galaxies.back().getMag() );
      //std::cout << " j in SourceMultiAnaGalaxy : " << j << std::endl;
			++j;
		}
	}

  std::cout << " closing file in SourceMultiAnaGalaxy : " << std::endl;
  color_cat.close();
	file_in.close();

	std::cout << galaxies.size() << " galaxies read in."<< std::endl;
	std::cout << "angular range x = " << rangex[0] << " to " << rangex[1] << "   y = " << rangey[0] << " to " << rangey[1]
  << std::endl;

	return;
}

void SourceMultiAnaGalaxy::assignParams(InputParams& params){
	if(!params.get("source_input_galaxy_file",input_gal_file)){
		  ERROR_MESSAGE();
		  std::cout << "parameter source_input_galaxy_file needs to be set in the parameter file "
				  << params.filename() << std::endl;
		  exit(0);
	  }
	if(!params.get("source_band",band)){
		std::cout << "ERROR: Must assign source_band in parameter file " << params.filename() << std::endl;
		std::cout << "Could be that specified band is not available " << std::endl;
    throw std::invalid_argument("band not supported");
	}
	if(!params.get("source_mag_limit",mag_limit)){
		std::cout << "ERROR: Must assign source_mag_limit in parameter file " << params.filename() << std::endl;
		exit(1);
	}

	if(!params.get("source_sb_limit",sb_limit))
		setSBlimit_magarcsec(30.);
	else
		sb_limit = pow(10,-0.4*(48.6+sb_limit))*pow(180*60*60/PI,2)/hplanck;
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
void SourceMultiAnaGalaxy::multiplier(
		PosType z                /// limiting redshift, only sources above this redshift are copied
		,PosType mag_cut         /// limiting magnitude, only sources with magnitudes below this limit will be copied
		,int multiplicity       /// the number of times each of these sources should be multiplied
		,Utilities::RandomNumbers_NR &ran    /// random number seed
		){
	unsigned long Nold = galaxies.size(),NtoAdd=0;
	PosType x1[2],x2[2],theta[2];

	for(unsigned long i=0;i<Nold;++i){
		if(galaxies[i].getTheta()[0] < x1[0]) x1[0] = galaxies[i].getTheta()[0];
		if(galaxies[i].getTheta()[0] > x2[0]) x2[0] = galaxies[i].getTheta()[0];
		if(galaxies[i].getTheta()[1] < x1[1]) x1[1] = galaxies[i].getTheta()[1];
		if(galaxies[i].getTheta()[1] > x2[1]) x2[1] = galaxies[i].getTheta()[1];

		if(galaxies[i].getZ() > z && galaxies[i].getMag() < mag_cut) ++NtoAdd;
	}

	NtoAdd = 0;
	for(unsigned long i=0;i<Nold;++i){
		if(galaxies[i].getZ() > z && galaxies[i].getMag() < mag_cut){

			for(int j=0;j<multiplicity;++j){
				theta[0] = x1[0] + (x2[0] - x1[0])*ran();
				theta[1] = x1[1] + (x2[1] - x1[1])*ran();

				galaxies.push_back(SourceOverzierPlus(galaxies[i].getMag(),galaxies[i].getMagBulge()
                ,galaxies[i].getReff(),galaxies[i].getRdisk(),ran()*PI,ran()*2*PI
					,Nold+NtoAdd,galaxies[i].getZ(),theta,ran));

				++NtoAdd;
			}
		}
	}

}

/// Sort the sources by redshift in assending order
void SourceMultiAnaGalaxy::sortInRedshift(){
  if(galaxies.size() < 2) return;
  std::sort(galaxies.begin(),galaxies.end(),redshiftcompare);
  delete searchtree;
  searchtree = new TreeSimpleVec<SourceOverzierPlus>(galaxies.data(),galaxies.size(),1,2,true,SourceOverzierPlus::getx);
}
// used in MultiSourceAnaGalaxy::sortInRedshift()
bool redshiftcompare(SourceOverzierPlus s1,SourceOverzierPlus s2){
    return (s1.getZ() < s2.getZ());
}
/// Sort the sources by magnitude in assending order
void SourceMultiAnaGalaxy::sortInMag(Band tmp_band){
  if(galaxies.size() < 2) return;
  
  /*if(tmp_band != band){
    for(size_t i = 0; i < galaxies.size() ; ++i) galaxies[i].setBand(tmp_band);
  }*/
  
  //std::sort(galaxies.begin(),galaxies.end(),magcompare);
  std::sort(galaxies.begin(),galaxies.end()
            ,[&tmp_band](SourceOverzierPlus s1,SourceOverzierPlus s2){
    return (s1.getMag(tmp_band) < s2.getMag(tmp_band));
  }
);
  
/*  if(tmp_band != band){
    for(size_t i = 0; i < galaxies.size() ; ++i) galaxies[i].setBand(band);
  }*/

  delete searchtree;
  searchtree = new TreeSimpleVec<SourceOverzierPlus>(galaxies.data(),galaxies.size(),1,2,true,SourceOverzierPlus::getx);
}
/// Sort the sources by magnitude in assending order
void SourceMultiAnaGalaxy::sortInID(){
  if(galaxies.size() < 2) return;
  std::sort(galaxies.begin(),galaxies.end(),
            [](SourceOverzierPlus s1,SourceOverzierPlus s2){return (s1.getID() < s2.getID());});

  delete searchtree;
  searchtree = new TreeSimpleVec<SourceOverzierPlus>(galaxies.data(),galaxies.size(),1,2,true,SourceOverzierPlus::getx);
}
// used in MultiSourceAnaGalaxy::sortInRedshift()
bool magcompare(SourceOverzierPlus s1,SourceOverzierPlus s2){
  return (s1.getMag() < s2.getMag());
}
// used in MultiSourceAnaGalaxy::sortInRedshift()
bool idcompare(SourceOverzierPlus s1,SourceOverzierPlus s2){
  return (s1.getID() < s2.getID());
}
/// Print info on current source parameters
void SourceMultiAnaGalaxy::printSource(){
	std::cout << "Overzier Galaxy Model" << std::endl;
	galaxies[index].printSource();
}
