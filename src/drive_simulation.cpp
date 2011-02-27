#include <stdlib.h>
#include <string>
#include <iostream>
#include <mesh_manager.h>
#include <messages.h>
#include <drive_simulation.h>

using namespace std;

int drive_simulation(Json::Value & config, Log & log, int argc, char * argv[] )
{
	start_message(log);

	// This information is hard coded for now, until we have an input reader:
	std::string meshFileName = config["mesh-database"].asString();
	std::string resultsFileName = config["results-database"].asString();
	const char* mesh_input_file_name = meshFileName.c_str();
	const char* mesh_output_file_name = resultsFileName.c_str();

	Mesh_Manager mesh_manager = Mesh_Manager(mesh_input_file_name, mesh_output_file_name);
	mesh_manager.read_mesh();

    stk::ParallelMachine parallel_machine = stk::parallel_machine_init(&argc, &argv);
    stk::mesh::STK_Mesh stk_mesh(parallel_machine, mesh_manager.dimension());

    mesh_manager.populate_STK_mesh(&stk_mesh);

    bool local_status = true ;
    local_status = mesh_manager.verify_coordinates_field(stk_mesh);
    stringstream oss;
    oss << "Verifying the STK mesh coordinates field: ";
    printStatus(log, local_status, &oss);

	//  CREATE BOGUS FIELDS FOR NOW the following is not linked to stk mesh stuff above yet

    const char * glo_var_1 = (const char *)"test_global_var_1";
    const char * glo_var_2 = (const char *)"test_global_var_2";
    mesh_manager.insert_global_var_name(glo_var_1);
    mesh_manager.insert_global_var_name(glo_var_2);
    const char * node_var_1 = (const char *)"test_nodal_var_1";
    const char * node_var_2 = (const char *)"test_nodal_var_2";
    mesh_manager.insert_nodal_var_name(node_var_1);
	mesh_manager.insert_nodal_var_name(node_var_2);
	const char * ele_var_1 = (const char *)"test_element_var_1";
	const char * ele_var_2 = (const char *)"test_element_var_2";
	mesh_manager.insert_element_var_name(ele_var_1);
	mesh_manager.insert_element_var_name(ele_var_2);

	const char * title = (const char*)"This is a test";
	mesh_manager.initialize_output(title);

	const int num_glo_vars = mesh_manager.num_global_variables();
	const int num_nod_vars = mesh_manager.num_nodal_variables();
	const int num_ele_vars = mesh_manager.num_element_variables();
	const int num_blocks = mesh_manager.num_blocks();

    const int num_nodes = mesh_manager.num_nodes();
    const int num_elem = mesh_manager.num_elem();
    const int CPU_word_size = mesh_manager.cpu_word_size();
    const int IO_word_size = mesh_manager.io_word_size();

	float glob_var_vals[num_glo_vars];
	float nodal_var_vals[num_nodes];

	int time_step = 1;
	int num_time_steps = 10;
	for (int i=0; i<num_time_steps; i++)
	{
		float time_value = (float)(i+1)/100.0;
		mesh_manager.write_time_step_info(time_step,time_value);
		for (int j=0; j<num_glo_vars; j++)
		{
			glob_var_vals[j] = (float)(j+1) * time_value;
		}
		mesh_manager.write_global_variables_to_output(time_step,time_value,glob_var_vals);
		for (int k=0; k<num_nod_vars; k++)
		{
			for (int j=0; j<num_nodes; j++)
			{
				nodal_var_vals[j] = (float)k + ((float)(j+1) * time_value);
			}
			mesh_manager.write_nodal_variable_to_output(time_step,time_value,nodal_var_vals,k+1);
		}
		for (int k=0; k<num_ele_vars; k++)
		{
			for (int j=0; j<num_blocks; j++)
			{
				const int num_elem_in_block = mesh_manager.num_elem_in_block(j);
				float elem_var_vals[num_elem_in_block];
				for (int m=0; m<num_elem_in_block; m++)
				{
					elem_var_vals[m] = (float)(k+1) + (float)(j+2) + ((float)(m+1)*time_value);
				}
				mesh_manager.write_element_variable_to_output(time_step,time_value,elem_var_vals,k+1,j);
			}
		}
		mesh_manager.update_output();
		time_step++;
	}

    mesh_manager.close_output_file();

	success_message(log);

	return 0;
}
