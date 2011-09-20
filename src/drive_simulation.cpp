#include <stdlib.h>
#include <string>
#include <iostream>
#include <mesh_manager.h>
#include <messages.h>
#include <drive_simulation.h>

using namespace std;

int drive_simulation(Json::Value & config, int argc, char * argv[])
{
  start_message();

  std::string meshFileName = config["mesh-database"].asString();
  std::string resultsFileName = config["results-database"].asString();
  const char* mesh_input_file_name = meshFileName.c_str();
  const char* mesh_output_file_name = resultsFileName.c_str();

  Mesh_Manager mesh_manager = Mesh_Manager(mesh_input_file_name,
    mesh_output_file_name);
  mesh_manager.read_mesh();

  // Note: we make this static so that drive_simulation can be called multiple times
  // FIXME: there's probably a better way to handle this.
  static stk::ParallelMachine parallel_machine = stk::parallel_machine_init(
    &argc, &argv);
  stk::mesh::STK_Mesh stk_mesh(parallel_machine, mesh_manager.dimension());

  //  CREATE BOGUS FIELDS FOR NOW
  stk::mesh::Part & universal = stk_mesh.my_fem_metaData.universal_part();

  stk::mesh::ScalarFieldType & temperature_field =
      stk_mesh.my_fem_metaData.declare_field<stk::mesh::ScalarFieldType>(
        "temperature");
  stk::mesh::put_field(temperature_field,
    stk_mesh.my_node_rank, universal);

  stk::mesh::VectorFieldType & velocity_field =
      stk_mesh.my_fem_metaData.declare_field<stk::mesh::VectorFieldType>(
        "velocity");
  stk::mesh::put_field(velocity_field,
    stk_mesh.my_node_rank, universal,
    stk_mesh.my_spatial_dimension);

  stk::mesh::ScalarFieldType & volume_field =
      stk_mesh.my_fem_metaData.declare_field<stk::mesh::ScalarFieldType>("volume"); //
  stk::mesh::put_field(volume_field,
    stk_mesh.my_elem_rank, universal);

  stk::mesh::VectorFieldType & acceleration_field =
      stk_mesh.my_fem_metaData.declare_field<stk::mesh::VectorFieldType>(
        "acceleration");
  stk::mesh::put_field(acceleration_field,
    stk_mesh.my_elem_rank, universal);

  mesh_manager.initialize_mesh_parts_and_commit(&stk_mesh);

  // This prints the field names that have been defined on all parts (universal)
  // It could be tailored to print specific parts like sidesets, etc.
  mesh_manager.print_field_info(&stk_mesh);
  // put the elements and connectivity in the stk_mesh
  mesh_manager.populate_mesh_elements(&stk_mesh);
  // put the coordinates as a field on the stk_mesh
  mesh_manager.populate_mesh_coordinates(&stk_mesh);
  bool local_status = true;
  local_status = mesh_manager.verify_coordinates_field(stk_mesh);
  stringstream oss;
  oss << "Verifying the STK mesh coordinates field: ";
  printStatus(local_status, &oss);

  // POPULATE BOGUS FIELDS FOR NOW
  // temperature
  mesh_manager.populate_bogus_scalar_field(&stk_mesh, temperature_field,
    stk_mesh.my_node_rank);
  // volume
  mesh_manager.populate_bogus_scalar_field(&stk_mesh, volume_field,
    stk_mesh.my_elem_rank);
  // velocity
  mesh_manager.populate_bogus_vector_field(&stk_mesh, velocity_field,
    stk_mesh.my_node_rank);
  // acceleration
  mesh_manager.populate_bogus_vector_field(&stk_mesh, acceleration_field,
    stk_mesh.my_elem_rank);

  const char * title = resultsFileName.c_str();
  mesh_manager.initialize_output(title, stk_mesh);

  int time_step = 1;
  int num_time_steps = 1;
  for (int i = 0; i < num_time_steps; i++)
  {
    float time_value = (float) (i + 1) / 100.0;
    mesh_manager.write_time_step_info(time_step, time_value);

    mesh_manager.write_nodal_scalar(time_step, time_value, stk_mesh,
      temperature_field);
    mesh_manager.write_nodal_vector(time_step, time_value, stk_mesh,
      stk_mesh.my_coordinates_field);
    mesh_manager.write_nodal_vector(time_step, time_value, stk_mesh,
      velocity_field);

    mesh_manager.write_element_scalar(time_step, time_value, stk_mesh,
      volume_field);
    mesh_manager.write_element_vector(time_step, time_value, stk_mesh,
      acceleration_field);

    mesh_manager.update_output();
    time_step++;
  }

  mesh_manager.close_output_file();

  success_message();
  return 0;
}
