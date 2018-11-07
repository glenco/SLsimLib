//
//  gadgetio.cpp
//  NR
//
//  Created by bmetcalf on 07.11.18.
//

#include "gadgetio.hpp"
#include <vector>
#include <stdio.h>
#include <iostream>


void load_snapshot_gadget2(char *fname
                           ,int nfiles
                           ,PosType **xp
                           ,std::vector<float> &masses
                           ,std::vector<float> &sizes
                           ,bool &multimass
                           )
{
  
  FILE *fd;
  char buf[200];
  int i, k, dummy, ntot_withmasses;
  int n, pc, pc_new, pc_sph;
  gadget_header header1;
  gadget_particle_data p;
  int Ngas = 0;
  int NumPart = 0;
  int id;
  
  std::vector<int> types;
  
  
  for(i = 0, pc = 1; i < nfiles; i++, pc = pc_new)
  {
    if(nfiles > 1)
      sprintf(buf, "%s.%d", fname, i);
    else
      sprintf(buf, "%s", fname);
    
    if(!(fd = fopen(buf, "r")))
    {
      printf("can't open file `%s`\n", buf);
      throw std::invalid_argument("no file");
    }
    
    printf("reading `%s' ...\n", buf);
    fflush(stdout);
    
    fread(&dummy, sizeof(dummy), 1, fd);
    fread(&header1, sizeof(header1), 1, fd);
    fread(&dummy, sizeof(dummy), 1, fd);
    
    if(nfiles == 1)
    {
      for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
        NumPart += header1.npart[k];
      Ngas = header1.npart[0];
    }
    else
    {
      for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
        NumPart += header1.npartTotal[k];
      Ngas = header1.npartTotal[0];
    }
    
    for(k = 0, ntot_withmasses = 0; k < 6; k++)
    {
      if(header1.mass[k] == 0)
        ntot_withmasses += header1.npart[k];
    }
    
    if(i == 0){
      // allocate memory
      try{
        //p.resize(NumPart);
        //ids.resize(NumPart);
        xp = Utilities::PosTypeMatrix(NumPart,3);
        masses.resize(NumPart);
        sizes.resize(NumPart,0);
        types.resize(NumPart);
      }catch (const std::bad_alloc& ba){
        std::cerr <<"ERROR: ";
        std::cerr << ba.what() << std::endl;
      }
    }
    
    float pos[3];
    float v_dummy[3];
    
    fread(&dummy, sizeof(dummy), 1, fd); // skip
    for(k = 0, pc_new = pc; k < 6; k++)
    {
      for(n = 0; n < header1.npart[k]; n++)
      {
        fread(pos, sizeof(float), 3, fd);
        xp[pc_new][0] = pos[0];
        xp[pc_new][1] = pos[0];
        xp[pc_new][2] = pos[0];
        pc_new++;
      }
    }
    fread(&dummy, sizeof(dummy), 1, fd); // skip
    
    fread(&dummy, sizeof(dummy), 1, fd); // skip
    for(k = 0, pc_new = pc; k < 6; k++)
    {
      for(n = 0; n < header1.npart[k]; n++)
      {
        fread(v_dummy, sizeof(float), 3, fd);
        pc_new++;
      }
    }
    fread(&dummy, sizeof(dummy), 1, fd); // skip
    
    
    fread(&dummy, sizeof(dummy), 1, fd); // skip
    for(k = 0, pc_new = pc; k < 6; k++)
    {
      for(n = 0; n < header1.npart[k]; n++)
      {
        fread(&id, sizeof(int), 1, fd);
        pc_new++;
      }
    }
    fread(&dummy, sizeof(dummy), 1, fd); // skip
    
    
    if(ntot_withmasses > 0)
      fread(&dummy, sizeof(dummy), 1, fd); // skip
    for(k = 0, pc_new = pc; k < 6; k++)
    {
      for(n = 0; n < header1.npart[k]; n++)
      {
        types[pc_new] = k;
        
        if(header1.mass[k] == 0)
          fread(&masses[pc_new], sizeof(float), 1, fd);
        else
          masses[pc_new] = header1.mass[k];
        pc_new++;
      }
    }
    
    if(ntot_withmasses > 0)
      fread(&dummy, sizeof(dummy), 1, fd); // skip
    
    
    float temp,rho,ne;
    
    if(header1.npart[0] > 0)
    {
      fread(&dummy, sizeof(dummy), 1, fd); // skip
      for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
      {
        fread(&temp, sizeof(float), 1, fd);
        pc_sph++;
      }
      fread(&dummy, sizeof(dummy), 1, fd); // skip
      
      fread(&dummy, sizeof(dummy), 1, fd); // skip
      for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
      {
        fread(&rho, sizeof(float), 1, fd);
        pc_sph++;
      }
      fread(&dummy, sizeof(dummy), 1, fd); // skip
      
      if(header1.flag_cooling)
      {
        fread(&dummy, sizeof(dummy), 1, fd); // skip
        for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
        {
          fread(&ne, sizeof(float), 1, fd);
          pc_sph++;
        }
        fread(&dummy, sizeof(dummy), 1, fd); // skip
      }
      else
        for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
        {
          ne = 1.0;
          pc_sph++;
        }
    }
    
    fclose(fd);
  }
}


/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
void load_snapshot_gadget2(char *fname, int nfiles, std::vector<gadget_particle_data> &p
                          ,std::vector<int> ids,gadget_header &header1
                          ,int &NumPart,int &Ngas
                          )
{
  FILE *fd;
  char buf[200];
  int i, k, dummy, ntot_withmasses;
  int n, pc, pc_new, pc_sph;

  
  for(i = 0, pc = 1; i < nfiles; i++, pc = pc_new)
  {
    if(nfiles > 1)
      sprintf(buf, "%s.%d", fname, i);
    else
      sprintf(buf, "%s", fname);
    
    if(!(fd = fopen(buf, "r")))
    {
      printf("can't open file `%s`\n", buf);
      throw std::invalid_argument("no file");
    }
    
    printf("reading `%s' ...\n", buf);
    fflush(stdout);
    
    fread(&dummy, sizeof(dummy), 1, fd);
    fread(&header1, sizeof(header1), 1, fd);
    fread(&dummy, sizeof(dummy), 1, fd);
    
    if(nfiles == 1)
    {
      for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
        NumPart += header1.npart[k];
      Ngas = header1.npart[0];
    }
    else
    {
      for(k = 0, NumPart = 0, ntot_withmasses = 0; k < 6; k++)
        NumPart += header1.npartTotal[k];
      Ngas = header1.npartTotal[0];
    }
    
    for(k = 0, ntot_withmasses = 0; k < 6; k++)
    {
      if(header1.mass[k] == 0)
        ntot_withmasses += header1.npart[k];
    }
    
    if(i == 0){
      // allocate memory
      try{
        p.resize(NumPart);
        ids.resize(NumPart);
      }catch (const std::bad_alloc& ba){
        std::cerr <<"ERROR: ";
        std::cerr << ba.what() << std::endl;
      }
    }
    
    
    gadget_particle_data *P = p.data() - 1;
    int *Id = ids.data() - 1;

    fread(&dummy, sizeof(dummy), 1, fd); // skip
    for(k = 0, pc_new = pc; k < 6; k++)
    {
      for(n = 0; n < header1.npart[k]; n++)
      {
        fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
        pc_new++;
      }
    }
    fread(&dummy, sizeof(dummy), 1, fd); // skip
    
    fread(&dummy, sizeof(dummy), 1, fd); // skip
    for(k = 0, pc_new = pc; k < 6; k++)
    {
      for(n = 0; n < header1.npart[k]; n++)
      {
        fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
        pc_new++;
      }
    }
    fread(&dummy, sizeof(dummy), 1, fd); // skip
    
    
    fread(&dummy, sizeof(dummy), 1, fd); // skip
    for(k = 0, pc_new = pc; k < 6; k++)
    {
      for(n = 0; n < header1.npart[k]; n++)
      {
        fread(&Id[pc_new], sizeof(int), 1, fd);
        pc_new++;
      }
    }
    fread(&dummy, sizeof(dummy), 1, fd); // skip
    
    
    if(ntot_withmasses > 0)
      fread(&dummy, sizeof(dummy), 1, fd); // skip
    for(k = 0, pc_new = pc; k < 6; k++)
    {
      for(n = 0; n < header1.npart[k]; n++)
      {
        P[pc_new].Type = k;
        
        if(header1.mass[k] == 0)
          fread(&P[pc_new].Mass, sizeof(float), 1, fd);
        else
          P[pc_new].Mass = header1.mass[k];
        pc_new++;
      }
    }
    
    if(ntot_withmasses > 0)
      fread(&dummy, sizeof(dummy), 1, fd); // skip
    
    
    if(header1.npart[0] > 0)
    {
      fread(&dummy, sizeof(dummy), 1, fd); // skip
      for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
      {
        fread(&P[pc_sph].U, sizeof(float), 1, fd);
        pc_sph++;
      }
      fread(&dummy, sizeof(dummy), 1, fd); // skip
      
      fread(&dummy, sizeof(dummy), 1, fd); // skip
      for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
      {
        fread(&P[pc_sph].Rho, sizeof(float), 1, fd);
        pc_sph++;
      }
      fread(&dummy, sizeof(dummy), 1, fd); // skip
      
      if(header1.flag_cooling)
      {
        fread(&dummy, sizeof(dummy), 1, fd); // skip
        for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
        {
          fread(&P[pc_sph].Ne, sizeof(float), 1, fd);
          pc_sph++;
        }
        fread(&dummy, sizeof(dummy), 1, fd); // skip
      }
      else
        for(n = 0, pc_sph = pc; n < header1.npart[0]; n++)
        {
          P[pc_sph].Ne = 1.0;
          pc_sph++;
        }
    }
    
    fclose(fd);
  }
}

