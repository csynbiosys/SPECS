/* Read parameter values from a file in format "<discarded> <value>"
 * 
 * M.Flores 2010
 */

#include "load_params.h"

void load_params(std::string filename, double* params, const int nParams) {
	//open file
	std::ifstream params_file(filename.data());
	if(!params_file) {
		std::cerr << "Failed to open the constants file " << filename.data() << std::endl;
		std::exit(EXIT_FAILURE);
	}

	const int buffer_length = 128;
	char buffer[buffer_length];
	char dummy[buffer_length];
	for(int i=0; i<nParams; ++i) {
		params_file.getline(buffer,buffer_length);
		std::sscanf(buffer,"%s %lf", dummy, &params[i]);
	}
	params_file.close();
	
	//return 0;
}
