#include <iostream>
#include <fstream>

#include <boost/program_options.hpp>

#include <json/json.h>

#include <drive_simulation.h>
#include <log.h>

namespace po = boost::program_options;

int main(int argc, char * argv[])
{
  // Declare the supported options.
  po::options_description desc("Allowed options");
  desc.add_options()("help,h", "produce help message")("input,i",
    po::value<std::string>(), "input file")("log,l", po::value<std::string>(), "output log file");

  // Parse the command line options
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  // Handle help requests
  if (vm.count("help"))
  {
    std::cerr << desc << "\n";
    return 1;
  }

  // Get the input file name
  std::string inputFileName = "input.json";
  if (vm.count("input"))
  {
    inputFileName = vm["input"].as<std::string>();
  }

  // Get the output log file name
  std::string logFileName = "output.log";
  if (vm.count("log"))
  {
    logFileName = vm["log"].as<std::string>();
  }

  // Setup the output log file. A file name of "-" causes the code
  // to use standard output.
  std::ofstream logFile;
  if (logFileName == "-")
  {
    Log::setLog(std::cout);
  }
  else
  {
    logFile.open(logFileName.c_str());
    Log::setLog(logFile);
  }
  // Set the error and warning streams too. We don't have options for these.
  Log::setWarning(std::cout);
  Log::setError(std::cerr);

  // Parse the input file
  log() << "  Using input file " << inputFileName << std::endl;
  std::ifstream inputFile;
  inputFile.open(inputFileName.c_str());
  Json::Value config;
  Json::Reader reader;
  const bool parsingSuccessful = reader.parse(inputFile, config);
  if (!parsingSuccessful)
  {
    std::cerr << "  Error parsing input file: " << inputFileName << std::endl;
    return 1;
  }

  // Run the actual simulation
  const int error_code = drive_simulation(config, argc, argv);
  return error_code;
}
