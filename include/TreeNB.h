 /*
 * Code Name:     treeNB.h
 * Programmer:    Ben Metcalf
 * Discription:                                                     
 * Comments:                           
 */

#ifndef treeNBtypes_declare
#define treeNBtypes_declare

#include <Tree.h>
#include <cosmo.h>

#ifndef pi
#define pi  3.141593
#endif

#ifndef treeNBdim
#define treeNBdim 2  // dimension of boxes in tree
#endif

/** type for particle positions and boundaries etc **/

#ifndef IndexType_declare
#define IndexType_declare
typedef unsigned long IndexType;
#endif

namespace Utilities{
	//PosType **PosTypeMatrix(long nrl, long nrh, long ncl, long nch);
	//void free_PosTypeMatrix(PosType **m, long nrl, long nrh, long ncl, long nch);
}
#endif
