#include <enums.h>
#include <map>

static bool string_maps_created = false;
static std::map<Physics_Type,std::string> physics_type_string;
static std::map<std::string,Physics_Type> string_physics_type;

void create_string_maps();
void create_string_maps()
{
  if (string_maps_created)
  {
    return;
  }
  string_maps_created = true;

  physics_type_string[PERIDYNAMICS]            = "PERIDYNAMICS";
  physics_type_string[FINITE_ELEMENT]          = "FINITE_ELEMENT";

  for (std::map<Physics_Type,std::string>::iterator pos = physics_type_string.begin(); pos != physics_type_string.end(); ++pos)
  {
    string_physics_type[pos->second] = pos->first;
  }
}

std::string tostring(const Physics_Type & physics_type)
{
  create_string_maps();
  std::map<Physics_Type,std::string>::iterator pos=physics_type_string.find(physics_type);
    if (pos == physics_type_string.end())
    {
      std::cerr << "Unknown physics type: " << physics_type << std::endl;
      exit(0);
    }
    return pos->second;
}

Physics_Type string_to_physics_type(const std::string & input_string)
{
  create_string_maps();
  std::map<std::string,Physics_Type,std::string >::iterator pos=string_physics_type.find(input_string);
  if (pos == string_physics_type.end())
  {
    std::cerr << "Unknown physics type: " << input_string << std::endl;
    std::cerr << "Valid options are: " << std::endl;
    for (std::map<std::string,Physics_Type,std::string >::iterator it = string_physics_type.begin(); it != string_physics_type.end(); ++it)
    {
      std::cerr << "  " << it->first << std::endl;
      exit(0);
    }
  }
  return pos->second;
}

