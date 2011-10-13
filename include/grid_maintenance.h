/*
 * grid_maintenance.h
 *
 *  Created on: Oct 12, 2011
 *      Author: bmetcalf
 */

GridHndl NewGrid(int Ngrid,double center[2],double range);
void FreeGrid(GridHndl grid);
void ReInitalizeGrid(GridHndl grid);
void TrimGrid(GridHndl grid,double highestres,Boolean useSB);
