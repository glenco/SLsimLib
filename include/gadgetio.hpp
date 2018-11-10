//
//  gadgetio.hpp
//  NR
//
//  Created by bmetcalf on 07.11.18.
//

#ifndef gadgetio_hpp
#define gadgetio_hpp

#include <vector>

struct gadget_header
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];  /* fills to 256 Bytes */
};

struct gadget_particle_data
{
  float Pos[3];
  float Vel[3];
  float Mass;
  int Type;
  
  float Rho, U, Temp, Ne;
};

void load_snapshot_gadget2(char *fname               /// base name of snapshot
                          ,int nfiles               /// number of files in snapshot
                          ,std::vector<gadget_particle_data> &p  /// particle information
                          ,std::vector<int> IDs     /// list of ID numbers in order
                          ,gadget_header &header1   /// header information from file
                          ,int &NumPart             /// total number of particles
                          ,int &Ngas                /// number of gas particles
                          );

void load_snapshot_gadget2(char *fname
                           ,int nfiles
                           ,std::vector<ParticleData<float> > &xp
                           ,std::vector<int> &ntypes
                           );
#endif /* gadgetio_hpp */
