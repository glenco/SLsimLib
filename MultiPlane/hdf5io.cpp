#include "particle_halo.h"

#include <stdexcept>

#ifdef ENABLE_HDF5
#include "H5Cpp.h"
#include <algorithm>
#include <array>
#include <string>
#include <vector>

bool MakeParticleLenses::readHDF5(){
  H5::Exception::dontPrint();
  H5::H5File file(filename.c_str(), H5F_ACC_RDONLY);
  H5::Group header = file.openGroup("Header");
  H5::Attribute redshift_attribute = header.openAttribute("Redshift");
  redshift_attribute.read(H5::PredType::NATIVE_DOUBLE, &z_original);

  std::array<double, 6> mass_table = {};
  const bool has_mass_table = header.attrExists("MassTable");
  if(has_mass_table){
    H5::Attribute mass_attribute = header.openAttribute("MassTable");
    mass_attribute.read(H5::PredType::NATIVE_DOUBLE, mass_table.data());
  }

  nparticles.assign(6, 0);
  masses.assign(6, 0);
  data.clear();
  const double length_scale = 1.0e-3 / (1.0 + z_original);

  for(int type = 0; type < 6; ++type){
    const std::string group_name = "PartType" + std::to_string(type);
    const htri_t group_exists = H5Lexists(file.getId(), group_name.c_str(), H5P_DEFAULT);
    if(group_exists < 0) throw std::runtime_error("Unable to inspect HDF5 particle groups.");
    if(group_exists == 0) continue;

    H5::Group particle_group = file.openGroup(group_name);
    H5::DataSet coordinates = particle_group.openDataSet("Coordinates");
    H5::DataSpace coordinate_space = coordinates.getSpace();
    hsize_t coordinate_dims[2] = {};
    if(coordinate_space.getSimpleExtentNdims() != 2
       || coordinate_space.getSimpleExtentDims(coordinate_dims) != 2
       || coordinate_dims[1] != 3){
      throw std::runtime_error("HDF5 Coordinates must have shape N x 3.");
    }

    const size_t count = static_cast<size_t>(coordinate_dims[0]);
    std::vector<double> position_values(count * 3);
    if(count > 0){
      coordinates.read(position_values.data(), H5::PredType::NATIVE_DOUBLE);
    }

    std::vector<double> mass_values(count);
    if(particle_group.nameExists("Masses")){
      H5::DataSet mass_dataset = particle_group.openDataSet("Masses");
      H5::DataSpace mass_space = mass_dataset.getSpace();
      hsize_t mass_dims[1] = {};
      if(mass_space.getSimpleExtentNdims() != 1
         || mass_space.getSimpleExtentDims(mass_dims) != 1
         || mass_dims[0] != coordinate_dims[0]){
        throw std::runtime_error("HDF5 Masses count does not match Coordinates.");
      }
      if(count > 0) mass_dataset.read(mass_values.data(), H5::PredType::NATIVE_DOUBLE);
    }else{
      if(!has_mass_table){
        throw std::runtime_error("HDF5 particle group has no Masses dataset or Header MassTable.");
      }
      std::fill(mass_values.begin(), mass_values.end(), mass_table[type]);
      masses[type] = static_cast<float>(mass_table[type]);
    }

    const size_t start = data.size();
    data.resize(start + count);
    for(size_t i = 0; i < count; ++i){
      ParticleType<float> &particle = data[start + i];
      for(int axis = 0; axis < 3; ++axis){
        particle[axis] = static_cast<float>(position_values[3 * i + axis] * length_scale);
      }
      particle.Mass = static_cast<float>(mass_values[i] * 1.0e10);
      particle.Size = 0;
      particle.type = type;
    }
    nparticles[type] = count;
  }

  if(data.empty()){
    throw std::runtime_error("HDF5 file contains no particles in PartType0-PartType5.");
  }

  size_t start = 0;
  for(int type = 0; type < 6; ++type){
    if(nparticles[type] > 0){
      LensHaloParticlesDep<ParticleType<float> >::calculate_smoothing(
        Nsmooth, data.data() + start, nparticles[type]);
    }
    start += nparticles[type];
  }
  return true;
}
#else
bool MakeParticleLenses::readHDF5(){
  throw std::runtime_error("HDF5 support is not enabled.");
}
#endif
