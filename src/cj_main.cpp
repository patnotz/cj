#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>

#include <json/json.h>

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

	// Get the input file name
	std::string inputFileName = "input.json";
	if ( vm.count("input")) {
		inputFileName = vm["input"].as<std::string>();
	}

	// Get the output log file name
	std::string logFileName = "output.log";
	if ( vm.count("log")) {
		logFileName = vm["log"].as<std::string>();
	}

	// Setup the output log file. A file name of "-" causes the code
	// to use standard output.
	if ( logFileName == "-") {
		Log::setLog(std::cout);
	} else {
		std::ofstream logFile;
		logFile.open(logFileName.c_str());
		Log::setLog(logFile);
	}
	// Set the error and warning streams too. We don't have options for these.
	Log::setWarning(std::cout);
	Log::setError(std::cerr);

	// Parse the input file
    log() << "Using input file " << inputFileName << std::endl;
	std::ifstream inputFile;
	inputFile.open(inputFileName.c_str());
	Json::Value config;
	Json::Reader reader;
	const bool parsingSuccessful = reader.parse(inputFile, config);
	if ( ! parsingSuccessful ) {
		std::cerr << "Error parsing input file: " << inputFileName << std::endl;
		return 1;
	}

	// If the user requested an alternate input mesh, override the input file:
	if ( vm.count("mesh")) {
		const std::string meshFileName = vm["mesh"].as<std::string>();
		config["mesh-database"] = meshFileName;
		log() << "Override input mesh file name: " << meshFileName << std::endl;
	}

	// If the user requested an alternate output mesh, override the input file:
	if ( vm.count("results")) {
		const std::string resultsFileName = vm["results"].as<std::string>();
		config["results-database"] = resultsFileName;
		log() << "Override output results file name: " << resultsFileName << std::endl;
	}

    // Run the actual simulation
	const int error_code = drive_simulation(config, argc, argv);
	return error_code;
}
