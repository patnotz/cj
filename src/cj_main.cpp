#include <stdlib.h>
#include <iostream>
#include <../include/mesh_manager.h>
#include <../include/messages.h>

//#include "exodusII.h"
//#include "netcdf.h"


using namespace std;

int main()
{
	start_message();

	// This information is hard coded for now, until we have an input reader:
	char* mesh_input_file_name = (char*)"mesh.g";
	char* mesh_output_file_name = (char*)"mesh.e";

	Mesh_Manager mesh_manager;
	mesh_manager.set_input_file(mesh_input_file_name);
	mesh_manager.read_mesh();
	mesh_manager.set_output_file(mesh_output_file_name);

	//  CREATE BOGUS OUTPUT FOR NOW

	char * glo_var_1 = (char *)"test_global_var_1";
	char * glo_var_2 = (char *)"test_global_var_2";
	mesh_manager.insert_global_var_name(glo_var_1);
	mesh_manager.insert_global_var_name(glo_var_2);
	char * node_var_1 = (char *)"test_nodal_var_1";
	char * node_var_2 = (char *)"test_nodal_var_2";
	mesh_manager.insert_nodal_var_name(node_var_1);
	mesh_manager.insert_nodal_var_name(node_var_2);
	char * ele_var_1 = (char *)"test_element_var_1";
	char * ele_var_2 = (char *)"test_element_var_2";
	mesh_manager.insert_element_var_name(ele_var_1);
	mesh_manager.insert_element_var_name(ele_var_2);

	char * title = (char*)"This is a test";
	mesh_manager.initialize_output(title);

	const int num_glo_vars = mesh_manager.num_global_variables();
	const int num_nod_vars = mesh_manager.num_nodal_variables();
	const int num_ele_vars = mesh_manager.num_element_variables();
	const int num_blocks = mesh_manager.num_blocks();

    const int num_nodes = mesh_manager.num_nodes();
    const int num_elem = mesh_manager.num_elem();
    const int CPU_word_size = mesh_manager.cpu_word_size();
    const int IO_word_size = mesh_manager.io_word_size();

	float * glob_var_vals = (float *) calloc (num_glo_vars, CPU_word_size);
	float * nodal_var_vals = (float *) calloc (num_nodes, CPU_word_size);

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
				float * elem_var_vals = (float *) calloc (num_elem_in_block, CPU_word_size);
				for (int m=0; m<num_elem_in_block; m++)
				{
					elem_var_vals[m] = (float)(k+1) + (float)(j+2) + ((float)(m+1)*time_value);
				}
				mesh_manager.write_element_variable_to_output(time_step,time_value,elem_var_vals,k+1,j);
				free(elem_var_vals);
			}
		}
//		mesh_manager.write_variables_to_output(time_step, time_value, glob_var_vals, nodal_var_vals, elem_var_vals);
		mesh_manager.update_output();
		time_step++;
	}

    mesh_manager.close_output_file();

    mesh_manager.gregs_output();

	success_message();
}


