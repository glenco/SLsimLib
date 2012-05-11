/*
 * tree_maintenance.h
 *
 *  Created on: Sep 29, 2011
 *      Author: bmetcalf
 */

#ifndef _tree_maintenance_declare_
#define _tree_maintenance_declare_

#include <Tree.h>

TreeHndl BuildTree(Point *xp,unsigned long Npoints,short my_median_cut = 1);
void _BuildTree(TreeHndl tree);
void FillTree(TreeHndl tree,Point *xp,unsigned long Npoints);
int AddPointsToTree(TreeHndl tree,Point *xpoint,unsigned long Nadd);
void _AddPoint(TreeHndl tree);
//unsigned long PruneTrees(TreeHndl i_tree,TreeHndl s_tree,double resolution,bool useSB);
unsigned long FreeBranchesBelow(TreeHndl i_tree,TreeHndl s_tree,KistHndl trashlist);
Point *RemoveLeafFromTree(TreeHndl tree,unsigned long *Npoints);

short emptyTree(TreeHndl tree);
short freeTree(TreeHndl tree);
void _freeBranches(TreeHndl tree,short child);
void _freeBranches_iter(TreeHndl tree);
void RebuildTreeFromList(TreeHndl tree);

#endif
