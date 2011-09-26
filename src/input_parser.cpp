#include <input_parser.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <analysis_model.h>
#include <Physics.h>
#include <messages.h>
#include <drive_simulation.h>
#include <boost/algorithm/string.hpp>
#include <enums.h>
#include <map>
#include <log.h>

using namespace std;

Input_Parser::Input_Parser(Json::Value & config, Simulation * simulation):
    my_config(config),
    my_simulation(simulation)
{
}

void Input_Parser::initialize()
{
  read_analysis_models();
  read_physics();
}

// Read in all the analysis model descriptions
void Input_Parser::read_analysis_models()
{
  Analysis_Model_Factory analysis_model_factory;
  std::map<std::string, Teuchos::RCP<Analysis_Model> > * analysis_models_map = my_simulation->get_analysis_models_map();
  ThrowRequire(analysis_models_map->empty());
  const Json::Value nullValue;
  unsigned int zero = 0;
  const int num_analysis_models = my_config["analysis-models"].size();
  for(unsigned int i = 0; i < num_analysis_models; ++i)
  {
    std::string analysis_model_name = my_config["analysis-models"].get(i,nullValue)["name"].asString();
    boost::to_upper(analysis_model_name);
    boost::algorithm::trim(analysis_model_name);
    std::string meshFileName = my_config["analysis-models"].get(i,nullValue)["mesh-database"].asString();
    std::string resultsFileName = my_config["analysis-models"].get(i,nullValue)["results-database"].asString();
    const char* mesh_input_file_name = meshFileName.c_str();
    const char* mesh_output_file_name = resultsFileName.c_str();
    Teuchos::RCP<Analysis_Model> analysis_model = analysis_model_factory.create(mesh_input_file_name, mesh_output_file_name,analysis_model_name);
    std::map<std::string, Teuchos::RCP<Analysis_Model> >::iterator pos=analysis_models_map->find(analysis_model_name);
    if (pos == analysis_models_map->end())
    {
      ThrowRequire(analysis_model_name == analysis_model->name()); //The string name should match its map index
      analysis_models_map->insert(std::pair<std::string, Teuchos::RCP<Analysis_Model> >(analysis_model_name,analysis_model));
    }
    else
    {
      std::cerr << "Input Parsing Error: analysis model names must be unique: " << analysis_model_name << std::endl;
      exit(1);
    }
  }
   ThrowRequire(!analysis_models_map->empty());
}

// Read in all the physics descriptions
void Input_Parser::read_physics()
{
  Physics_Factory physics_factory;
  const Json::Value nullValue;
  std::map<std::string, Teuchos::RCP<Analysis_Model> > * analysis_models_map = my_simulation->get_analysis_models_map();
  std::map<std::string, Teuchos::RCP<Physics> > * physics_map = my_simulation->get_physics_map();
  ThrowRequire(!analysis_models_map->empty());
  const int num_physics = my_config["physics"].size();
  for(unsigned int i = 0; i < num_physics; ++i)
  {
    std::string physics_name = my_config["physics"].get(i,nullValue)["name"].asString();
    boost::to_upper(physics_name);
    boost::algorithm::trim(physics_name);
    std::string physics_type_str = my_config["physics"].get(i,nullValue)["type"].asString();
    boost::to_upper(physics_type_str);
    const Physics_Type physics_type = string_to_physics_type(physics_type_str);
    std::string analysis_model_name = my_config["physics"].get(i,nullValue)["analysis-model"].asString();
    if(analysis_model_name=="")
    {
      std::cerr << "Input Parsing Error: every physics must have an analysis_model: " << physics_name << std::endl;
      exit(1);
    }
    boost::to_upper(analysis_model_name);
    boost::algorithm::trim(analysis_model_name);
    std::map<std::string, Teuchos::RCP<Analysis_Model> >::iterator pos=analysis_models_map->find(analysis_model_name);
    if (pos == analysis_models_map->end())
    {
      std::cerr << "Input Parsing Error: analysis model does not exist: " << analysis_model_name << std::endl;
      exit(1);
    }
    else
    {
      Teuchos::RCP<Physics> physics = physics_factory.create(analysis_model_name,physics_type,physics_name);
      std::map<std::string, Teuchos::RCP<Physics> >::iterator pos=physics_map->find(physics_name);
      if (pos == physics_map->end())
      {
        ThrowRequire(physics_name==physics->name()); // map index should match the physics name
        physics_map->insert(std::pair<std::string, Teuchos::RCP<Physics> >(physics_name,physics));
      }
      else
      {
        std::cerr << "Input Parsing Error: physics names must be unique: " << physics_name << std::endl;
        exit(1);
      }
    }
  }
  ThrowRequire(!physics_map->empty());

}
