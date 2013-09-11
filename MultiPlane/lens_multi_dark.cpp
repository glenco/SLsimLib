#include "../include/lens.h"
#include "../include/MOKAlens.h"

#include <fstream>
#include <stdexcept>

namespace
{
	const char* char_comment = "#";
	const char* char_whitespace = " \t";
}

void Lens::readMultiDark()
{
	std::ifstream list(main_input_file.c_str());
	if(!list.good())
		throw new std::runtime_error("Could not open MultiDark list file" + main_input_file + "!");
	
	std::cout << "MultiDark files:" << std::endl;
	
	std::string line;
	while(std::getline(list, line))
	{
		// skip empty lines
		if(line.empty())
			continue;
		
		std::ifstream::pos_type begin, end, comment;
		
		comment = line.find_first_of(char_comment);
		begin = line.find_first_not_of(char_whitespace);
		end = line.find_last_not_of(char_whitespace, comment - (std::ifstream::pos_type)1);
		
		// skip empty filenames
		if(begin >= end)
			continue;
		
		std::string mokafile = line.substr(begin, 1 + end - begin);
		
		std::cout << "- " << mokafile << std::endl;
		
		// create the MOKA halo
		main_halos.push_back(new LensHaloMOKA(mokafile));
	}
}
