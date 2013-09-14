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
		main_halos.push_back(new LensHaloMOKA(mokafile));
	}
}
