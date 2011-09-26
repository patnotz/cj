#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <map>
#include <analysis_model.h>
#include <physics.h>

class Simulation
{
public:
  static Simulation * Instance();

  std::map< std::string, Teuchos::RCP<Analysis_Model> > * get_analysis_models_map(){return &my_analysis_models_map;}
  std::map< std::string, Teuchos::RCP<Physics> > * get_physics_map(){return &my_physics_map;}
  void print_analysis_models();
  void print_physics();
  void read_meshes();
  void create_physics_stk_meshes(int argc, char * argv[]);
  void initialize_fields();

  void create_bogus_fields_and_output(int argc, char * argv[]);

private:
  Simulation();
  std::map< std::string, Teuchos::RCP<Analysis_Model> > my_analysis_models_map;
  std::map< std::string, Teuchos::RCP<Physics> > my_physics_map;
  static Simulation * my_simulation_ptr;
};


#endif /* SIMULATION_H_ */
