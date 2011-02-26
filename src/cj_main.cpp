#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>

#include <drive_simulation.h>
#include <log.h>

namespace po = boost::program_options;

int main( int argc, char * argv[] )
{
	// Declare the supported options.
	po::options_description desc("Allowed options");
	desc.add_options()
	    ("help,h", "produce help message")
	    ("input,i", po::value<std::string>(), "input file")
	    ("mesh,m", po::value<std::string>(), "input mesh database (overrides input file)")
	    ("results,r", po::value<std::string>(), "output results database (overrides input file)")
	    ("log,l", po::value<std::string>(), "output log file")
	;

	// Parse the command line options
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	// Handle help requests
	if (vm.count("help")) {
	    std::cerr << desc << "\n";
	    return 1;
	}

	// Setup the output log file. No log option or a  file name of "-"
	// causes the code to use standard output.
	Log *logPtr = NULL;
	std::ofstream logFile;
	if ( (! vm.count("log")) || vm["log"].as<std::string>() == "-") {
		logPtr = & std::cout;
	} else {
		std::string logFileName = vm["log"].as<std::string>();
		logFile.open(logFileName.c_str());
		logPtr = & logFile;
	}
	Log & log = *logPtr;

	// Read the input file
	std::ifstream inputFile;
	if ( ! vm.count("input")) {
	    std::cerr << "Input file not set.\n";
	    // return 1; // this will become and error soon
	}
	const std::string inputFileName = vm["input"].as<std::string>();
    log << "Using input file " << inputFileName << ".\n";
	inputFile.open(inputFileName.c_str());

	// Run the actual simulation
	const int error_code = drive_simulation(log, argc, argv);
	return error_code;
}
