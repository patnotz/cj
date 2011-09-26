#include <stdlib.h>
#include <string>
#include <iostream>
#include <analysis_model.h>
#include <Physics.h>
#include <messages.h>
#include <drive_simulation.h>
#include <simulation.h>
#include <input_parser.h>
#include <boost/algorithm/string.hpp>
#include <enums.h>
#include <map>
#include <log.h>

using namespace std;

static std::map< std::string, Teuchos::RCP<Physics> > physics_map;

int drive_simulation(Json::Value & config, int argc, char * argv[])
{
  start_message();

  // Create a simulation that will manage the meshes, equations, physics, etc
  Simulation * simulation = Simulation::Instance(); // only one simulation allowed
  Input_Parser input_parser(config, simulation);
  input_parser.initialize();

  // Now need to read in all the meshes in the simulation
  simulation->read_meshes();
  simulation->create_physics_stk_meshes(argc,argv);
  simulation->initialize_fields();

  success_message();
  return 0;
}
