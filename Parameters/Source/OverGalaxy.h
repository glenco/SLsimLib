#ifndef PARAMETERS_OVERGALAXY_H_
#define PARAMETERS_OVERGALAXY_H_

#include "../../include/parameters.h"

#include "../../Galaxies/overzier.h"

template<>
struct Parameters<OverGalaxy>
{
	/// haloID
	unsigned long haloID;
	
	/// redshift
	double z;
	
	/// total magnitude
	double mag;
	
	/// bulge to total ratio
	double BtoT;
	
	/// bulge half light radius (arcs)
	double Reff;
	
	/// disk scale height (arcs)
	double Rh;
	
	/// position angle (radians)
	double PA;
	
	/// inclination of disk (radians)
	double inclination;
	
	/// position on the sky
	double theta[2];	
};

void operator<<(Parameters<OverGalaxy>& p, const OverGalaxy& g)
{
	p.haloID = g.haloID;
	p.z = g.z;
	p.mag = g.getMag();
	p.BtoT = g.getBtoT();
	p.Reff = g.getReff();
	p.Rh = g.getRh();
	p.PA = g.getPA();
	p.inclination = g.getInclination();
	p.theta[0] = g.theta[0];
	p.theta[1] = g.theta[1];
}

void operator>>(const Parameters<OverGalaxy>& p, OverGalaxy& g)
{
	g.setInternals(p.mag, p.BtoT, p.Reff, p.Rh, p.PA, p.inclination, p.haloID, p.z, p.theta);
}

#endif
