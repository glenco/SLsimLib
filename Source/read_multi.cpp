#include "../include/multi_source.h"
#include "../include/overzier_source.h"

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

void MultiSource::readGalaxyFile(std::string filename, Band band, double mag_limit)
{
	std::size_t count = 0;
	
	//char c='0';
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
	
	std::ifstream file_in(filename.c_str());
	if(!file_in){
		std::cerr << "Can't open file " << filename << std::endl;
		exit(1);
	}
	
	std::cout << "Reading from galaxy data file " << filename.c_str() << std::endl;
	//file_in >> Ngalaxies;
	//std::cout << "Number of source galaxies: " << Ngalaxies << std::endl;
	
	i=0;
	while(file_in.peek() == '#'){
		file_in.ignore(10000,'\n');
		++i;
	}
	std::cout << "skipped "<< i << " comment lines in " << filename << std::endl;
	
	double theta[2] = {0.0,0.0};
	
	const int ncolumns = 31;
	
	void* addr[ncolumns];
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
	
	double mag,mag_bulge;
	
	// read in data
	for(i=0,j=0 ; ; ++i){
		myline.clear();
		getline(file_in,myline);
		
		if(myline[0] == '#')
			break;
		
		for(int l=0;l<ncolumns; l++){
			std::size_t pos = myline.find(f);
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
				std::cerr << "Requested band is not an available option." << std::endl;
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
			
			/***************************/
			OverzierSource* over = new OverzierSource(mag, pow(10, -(mag_bulge-mag)/2.5), Ref, Rh, pa, inclination, HaloID, z_cosm, theta);
			
			over->setUMag(SDSS_u);
			over->setGMag(SDSS_g);
			over->setRMag(SDSS_r);
			over->setIMag(SDSS_i);
			over->setZMag(SDSS_z);
			over->setJMag(J_band);
			over->setHMag(H_band);
			over->setKMag(Ks_band);
			
			// add galaxy to list of sources
			addInternal(over, typeid(OverzierSource), true);
			
			// increase counter
			++count;
			
			//std::cout << "z:" << z_cosm << " mag " << SDSS_u << " Bulge to total " << pow(10,-(SDSS_u_Bulge-SDSS_u)/2.5)
			//		<< " bulge size arcsec " << Ref  << " disk size arcsec " << pa << " position angle " << pa << " inclination " << inclination
			//		<< " theta = " << theta[0] << " " << theta[1] << std::endl;
			
			++j;
		}
	}
	
	file_in.close();
	
	std::cout << "Done reading " << count << " galaxies." << std::endl;
	
	return;
}
