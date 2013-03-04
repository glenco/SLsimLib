/*
 * tree_maintenance.h
 *
 *  Created on: Sep 29, 2011
 *      Author: bmetcalf
 */

#ifndef _tree_maintenance_declare_
#define _tree_maintenance_declare_

#include "Tree.h"

//unsigned long PruneTrees(TreeHndl i_tree,TreeHndl s_tree,double resolution,bool useSB);
unsigned long FreeBranchesBelow(TreeHndl i_tree,TreeHndl s_tree,KistHndl trashlist);
Point *RemoveLeafFromTree(TreeHndl tree,unsigned long *Npoints);

#endif
