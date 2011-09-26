#include <simulation.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <analysis_model.h>
#include <messages.h>
#include <drive_simulation.h>
#include <boost/algorithm/string.hpp>
#include <enums.h>
#include <map>
#include <log.h>
#include <physics.h>
#include <Teuchos_RCP.hpp>

using namespace std;

Simulation * Simulation::my_simulation_ptr = 0;

Simulation * Simulation::Instance()
{
  if (!my_simulation_ptr)   // Only allow one instance of class to be generated.
        my_simulation_ptr = new Simulation;
     return my_simulation_ptr;
}


Simulation::Simulation()
{
  my_analysis_models_map.clear();
  my_physics_map.clear();
}

void Simulation::print_analysis_models()
{
  log() << endl << "Analysis models: " << endl << endl;
  for(std::map<std::string, Teuchos::RCP<Analysis_Model> >::iterator pos=my_analysis_models_map.begin(); pos!=my_analysis_models_map.end();++pos)
  {
    log() << "Name: " << pos->second->name() << endl
        << "    Input file: " << pos->second->input_file_name() << endl
        << "    Output file: " << pos->second->output_file_name() << endl;
  }

}
void Simulation::print_physics()
{
  log() << endl << "Phyics: " << endl << endl;
  for(std::map<std::string, Teuchos::RCP<Physics> >::iterator pos=my_physics_map.begin(); pos!=my_physics_map.end();++pos)
  {
    log() << "Name: " << pos->second->name() << endl
          << "    Physics type: " << tostring(pos->second->physics_type())  << endl
          << "    Analysis model name: " << pos->second->analysis_model_str() << endl;
  }
}

void Simulation::read_meshes()
{
  std::map<std::string, Teuchos::RCP<Analysis_Model> >::iterator pos=my_analysis_models_map.begin();
  for(; pos!=my_analysis_models_map.end();++pos)
  {
     pos->second->read_mesh();
  }
}

void Simulation::create_physics_stk_meshes(int argc, char * argv[])
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Simulation::create_physics_stk_meshes()";
  oss << "Creating the parallel machine and the stk mesh";
  progress_message(&oss, method_name);
#endif

  // create an stk mesh instance and assign the physics stk_mesh pointer to it
  static stk::ParallelMachine parallel_machine = stk::parallel_machine_init(
        &argc, &argv);

  std::map<std::string, Teuchos::RCP<Physics> >::iterator pos=my_physics_map.begin();
  for(; pos!=my_physics_map.end();++pos)
  {
    const std::string analysis_model_str = pos->second->analysis_model_str();
    //get a pointer to the analysis_model for this physics
    std::map<std::string, Teuchos::RCP<Analysis_Model> >::iterator pos_a=my_analysis_models_map.find(analysis_model_str);
    if(pos_a==my_analysis_models_map.end())
    {
      cerr << "Error: could not find this physics analysis model in the simulation's map: physics name = " << pos->first << " analysis model = " << analysis_model_str << endl;
      exit(1);
    }
    const int dimension = pos_a->second->dimension();
#ifdef DEBUG_OUTPUT
  oss << "Creating an stk mesh for physics: " << pos->second->name();
  sub_progress_message(&oss);
#endif
    pos->second->stk_mesh_ptr() = Teuchos::rcp(new stk::mesh::STK_Mesh(parallel_machine, dimension));
  }
}
void Simulation::initialize_fields()
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Simulation::initialize_fields()";
  oss << "Initializing fields for each physics";
  progress_message(&oss, method_name);
#endif
  // create the fields for each physics
  std::map<std::string, Teuchos::RCP<Physics> >::iterator pos=my_physics_map.begin();
  for(; pos!=my_physics_map.end();++pos)
  {
    const std::string analysis_model_str = pos->second->analysis_model_str();
    //get a pointer to the analysis_model for this physics
    std::map<std::string, Teuchos::RCP<Analysis_Model> >::iterator pos_a=my_analysis_models_map.find(analysis_model_str);
    if(pos_a==my_analysis_models_map.end())
    {
      cerr << "Error: could not find this physics analysis model in the simulation's map: physics name = " << pos->first << " analysis model = " << analysis_model_str << endl;
      exit(1);
    }
    pos->second->initialize_fields(*pos_a->second);
  }
}
