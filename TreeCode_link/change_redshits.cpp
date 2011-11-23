/*
 * chang_redshits.c
 *
 *  Created on: Jun 15, 2010
 *      Author: R.B. Metcalf
 *
 *      changes redshift of source and/or lens without having to recalculate
 *      deflection angles
 *
 *      the image point remain at the same physical points on the lens plane
 *      the physical size of the source remains constant
 */

/*#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <Tree.h>
#include <tree_maintenance.h>
#include <analytic_lens.h>
*/
//#include "../TreeCode/TreeNB.h"

#include <slsimlib.h>

extern COSMOLOGY cosmo;
/** \ingroup ChangeLens
 *
 * \brief SHOULD BE TESTED Changes source and/or source redshits and recalculates the
 * the source points.  DOES NOT YET RECALCULATE SOURCE TREE
 */
void change_redshifts(TreeHndl i_tree,TreeHndl s_tree,AnaLens *lens,double z_source
		,double z_lens){
	double oldSigma=0,factor=0;

	oldSigma=lens->Sigma_crit;
	lens->source_r *= cosmo.angDist(0,z_source)/cosmo.angDist(0,lens->zsource);
	lens->zlens=z_lens;
	lens->zsource=z_source;

	lens->Sigma_crit = cosmo.angDist(0,lens->zsource)
			/cosmo.angDist(lens->zlens,lens->zsource)/cosmo.angDist(0,lens->zlens)/4/pi/Grav;

	factor=oldSigma/lens->Sigma_crit;

	lens->MpcToAsec=60*60*180*(1+lens->zsource)/pi/cosmo.angDist(0,lens->zlens);
	lens->host_ro=4*pi*pow(lens->host_sigma/2.99792e5,2)*cosmo.angDist(0,lens->zlens)
			*cosmo.angDist(lens->zlens,lens->zsource)
			  /cosmo.angDist(0,lens->zsource)/(1+lens->zlens);

	if(lens->Sigma_crit == oldSigma) return ;

	// shift the source points
	MoveToTopList(i_tree->pointlist);
	do{
		i_tree->pointlist->current->image->x[0] = (1-factor)*i_tree->pointlist->current->x[0]
		                                     + factor*i_tree->pointlist->current->image->x[0];
		i_tree->pointlist->current->image->x[1] = (1-factor)*i_tree->pointlist->current->x[1]
		                                     + factor*i_tree->pointlist->current->image->x[1];
	}while(MoveDownList(i_tree->pointlist));

	// rebuild tree on source plane
	moveTop(s_tree);
	_freeBranches(s_tree,0);
	assert(s_tree->Nbranches==0);
	assert(s_tree->pointlist->Npoints > 0);
	MoveToTopList(s_tree->pointlist);
	_BuildTree(s_tree);

	return ;
}