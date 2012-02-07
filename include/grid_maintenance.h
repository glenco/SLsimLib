/*
 * grid_maintenance.h
 *
 *  Created on: Oct 12, 2011
 *      Author: bmetcalf
 */

#include <model.h>
#include <lens.h>

GridHndl NewGrid(LensHndl lens,int Ngrid,double center[2],double range);
void FreeGrid(GridHndl grid);
void ReInitializeGrid(LensHndl lens,GridHndl grid);
void TrimGrid(GridHndl grid,double highestres,bool useSB);
void RefreshSurfaceBrightnesses(SourceHndl source,GridHndl grid);
unsigned long NumberOfPoints(GridHndl grid);
