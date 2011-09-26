#ifndef INPUT_PARSER_H_
#define INPUT_PARSER_H_

#include <stdlib.h>
#include <string>
#include <iostream>
#include <analysis_model.h>
#include <Physics.h>
#include <boost/algorithm/string.hpp>
#include <enums.h>
#include <map>
#include <log.h>
#include <simulation.h>
#include <json/json.h>

class Input_Parser
{
public:
  Input_Parser(Json::Value & config, Simulation * simulation);
  ~Input_Parser(){};

  void initialize();
  void read_analysis_models();
  void read_physics();

private:

  Json::Value my_config;
  Simulation * my_simulation;
};


#endif /* INPUT_PARSER_H_ */
