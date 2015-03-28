#include "../include/lens.h"
#include "../include/MOKAlens.h"

#include <fstream>
#include <stdexcept>

namespace
{
	const char* char_comment = "#";
	const char* char_whitespace = " \t";
}

void Lens::readPixelizedDensity()
{
	std::ifstream list(pixel_map_input_file.c_str());
	if(!list.good())
	{
		std::cerr << "Cannot open file " << pixel_map_input_file << std::endl;
		throw std::runtime_error("Could not open PixelDMap list file" + pixel_map_input_file + "!");
	}
	
  if(pixel_map_input_file.find(".fits") != pixel_map_input_file.npos){
    std::cout << "inputing PixelDensityMap file: " << pixel_map_input_file << std::endl;
    main_halos.push_back(new LensHaloMassMap(pixel_map_input_file, pix_map, cosmo));
  }else{
    
    std::cout << "reading PixelDMap files: " << pixel_map_input_file << std::endl;
    
    std::string line;
    while(std::getline(list, line))
    {
      // skip empty lines
      if(line.empty())
        continue;
      
      std::string::size_type comment = line.find_first_of(char_comment);
      
      if(comment == 0)
        continue;
      if(comment != std::string::npos)
        line = line.substr(0, comment - 1);
      
      std::string::size_type begin = line.find_first_not_of(char_whitespace);
      std::string::size_type end = line.find_last_not_of(char_whitespace);
      
      // skip empty filenames
      if(end <= begin)
        continue;
      
      std::string mokafile = line.substr(begin, end - begin + 1);
      
      std::cout << "- " << mokafile << std::endl;
      
      // create the MOKA halo
      main_halos.push_back(new LensHaloMassMap(mokafile, pix_map, cosmo));
    }
  }
}
