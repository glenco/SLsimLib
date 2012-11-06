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

#include "slsimlib.h"

/** \ingroup ChangeLens
 *
 * \brief SHOULD BE TESTED Changes source and/or source redshits and recalculates the
 * the source points.  DOES NOT YET RECALCULATE SOURCE TREE
 */
/*
void Model::change_redshifts(TreeHndl i_tree,TreeHndl s_tree,double z_source
		,double z_lens){
	double oldSigma=0,factor=0;

	oldSigma=lens->Sigma_crit;
	source->source_r *= cosmo->angDist(0,z_source)/cosmo->angDist(0,source->zsource);

	// chnage the redshifts
	lens->zlens=z_lens;
	source->zsource=z_source;

	setInternal();

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
}*/
