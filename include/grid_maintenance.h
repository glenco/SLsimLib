/*
 * grid_maintenance.h
 *
 *  Created on: Oct 12, 2011
 *      Author: bmetcalf
 */

GridHndl NewGrid(int Ngrid,double center[2],double range);
void FreeGrid(GridHndl grid);
void ReInitializeGrid(GridHndl grid);
void TrimGrid(GridHndl grid,double highestres,bool useSB);
void RefreshSurfaceBrightnesses(GridHndl grid,AnaLens *lens);
unsigned long NumberOfPoints(GridHndl grid);
