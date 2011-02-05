/*
 * mesh_manager.cpp
 *
 *  Created on: Jan 30, 2011
 *      Author: dzturne1
 */

#include <iostream>
#include <string>
#include <sstream>
#include <time.h>

#include <../include/mesh_manager.h>
#include <../include/messages.h>
#include <../include/main.h>

#include "exodusII.h"
#include "netcdf.h"

using namespace std;

void
Mesh_Manager::read_mesh()
{
	stringstream oss;
	int error;
	float version;

	my_CPU_word_size = 0;
	my_IO_word_size = 0;
	/*open exodus II file */
	my_input_exoid = ex_open(my_input_file_name,EX_READ,&my_CPU_word_size,&my_IO_word_size,&version);
	if(my_input_exoid<0)
	{
		oss << "Reading mesh failure.";
		sub_progress_message(&oss);
		exit(1);
	}

	initialize_read();
	import_nodes();
	import_elem_map();
	import_blocks();
    import_connectivities();
    import_node_sets();
    import_side_sets();

    error = ex_close(my_input_exoid);
}

void
Mesh_Manager::initialize_read()
{
	  string method_name = "Mesh_Manager::initialize_read()";
	  int error;
	  stringstream oss;

#ifdef DEBUG_OUTPUT
	  oss << "Reading file: " << my_input_file_name;
	  progress_message(&oss, method_name);
#endif
	  /* read database parameters */
	  char title[MAX_LINE_LENGTH+1];

	  error = ex_get_init (my_input_exoid, title, &my_num_dim, &my_num_nodes, &my_num_elem,
	  &my_num_elem_blk, &my_num_node_sets, &my_num_side_sets);

	  cout << endl <<
			  "  ----------------------------------------------------------------------------" << endl <<
			  "    Title: " << title << endl <<
			  "    Spatial dimension: " << my_num_dim << endl <<
			  "    " << my_num_nodes << " node(s)   " << my_num_elem << " element(s)   " << my_num_elem_blk
			  << " block(s)   " << my_num_node_sets << " node set(s)   " <<  my_num_side_sets << " side set(s)" << endl <<
			  "  ----------------------------------------------------------------------------" << endl;
}

void
Mesh_Manager::import_nodes()
{
	int error;
	/* read nodal coordinates values and names from database */
	my_x = (float *) calloc(my_num_nodes, sizeof(float));
	my_y = (float *) calloc(my_num_nodes, sizeof(float));
	if (my_num_dim >= 3)
		my_z = (float *) calloc(my_num_nodes, sizeof(float));
	else
		my_z = 0;
	error = ex_get_coord (my_input_exoid, my_x, my_y, my_z);

//  free (x);
//  free (y);
//  if (my_num_dim >= 3)
//  free (z);
}

void
Mesh_Manager::import_elem_map()
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::import_elem_map()";
	stringstream oss;
	oss << "Reading element map";
	progress_message(&oss,method_name);
#endif

	int error;
	/* read element order map */
	my_elem_map = (int *) calloc(my_num_elem, sizeof(int));
	error = ex_get_map (my_input_exoid, my_elem_map);
}

void
Mesh_Manager::import_blocks()
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::import_blocks()";
	stringstream oss;
	oss << "Reading blocks";
	progress_message(&oss,method_name);
#endif

	char elem_type[MAX_STR_LENGTH+1];
	int error;

	/* read element block parameters */
	my_block_ids = (int *) calloc(my_num_elem_blk, sizeof(int));
	my_num_elem_in_block = (int *) calloc(my_num_elem_blk, sizeof(int));
	my_num_nodes_per_elem = (int *) calloc(my_num_elem_blk, sizeof(int));
	my_num_attr = (int *) calloc(my_num_elem_blk, sizeof(int));
	error = ex_get_elem_blk_ids (my_input_exoid, my_block_ids);
	for (int i=0; i<my_num_elem_blk; i++)
	{
		error = ex_get_elem_block (my_input_exoid, my_block_ids[i], elem_type,
				&(my_num_elem_in_block[i]),
				&(my_num_nodes_per_elem[i]), &(my_num_attr[i]));

        my_elem_types.insert(pair<int,string>(my_block_ids[i],elem_type));

#ifdef DEBUG_OUTPUT
		oss << "Element block: " << my_block_ids[i];
		sub_progress_message(&oss);
		oss << "Element type: " << elem_type;
		sub_sub_progress_message(&oss);
		oss << "Number of elements in block: " << my_num_elem_in_block[i];
		sub_sub_progress_message(&oss);
		oss << "Number of nodes per element: " << my_num_nodes_per_elem[i];
		sub_sub_progress_message(&oss);
		oss << "Number of attributes: " << my_num_attr[i];
		sub_sub_progress_message(&oss);
#endif
	};
}

void
Mesh_Manager::import_connectivities()
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::import_connectivity()";
	stringstream oss;
#endif

	int error;
	int * connectivity;

	/* read element connectivity */
	for (int i=0; i<my_num_elem_blk; i++)
	{
#ifdef DEBUG_OUTPUT
		oss << "Reading connectivity for block " << my_block_ids[i];
		progress_message(&oss, method_name);
#endif
		connectivity = (int *) calloc((my_num_nodes_per_elem[i] * my_num_elem_in_block[i]),
				sizeof(int));
		error = ex_get_elem_conn (my_input_exoid, my_block_ids[i], connectivity);
        my_connectivities.insert(pair<int,int*>(my_block_ids[i],connectivity));
		// print_connectivity(my_block_ids[i]);
	}
}

void
Mesh_Manager::print_connectivity(int block_id)
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::print_connectivity()";
	stringstream oss;
	oss << "Block ID: " << block_id;
	progress_message(&oss,method_name);
#endif
	if(my_connectivities.find(block_id)==my_connectivities.end())
		return;

	int * connectivity = my_connectivities.find(block_id)->second;

    int num_elem = my_num_elem_in_block[block_id - 1];
	int num_nodes_per_elem = my_num_nodes_per_elem[block_id - 1];
	for(int i = 0;i < num_elem;++i)
	{
		cout << "Element " << i << ": ";
		for(int j = 0;j < num_nodes_per_elem;++j)
		{
			const int index = i*num_nodes_per_elem + j;
			cout << connectivity[index] << " ";
		}
		cout << endl;
	}
}

void
Mesh_Manager::import_node_sets()
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::import_node_sets()";
	stringstream oss;
	oss << "Reading node sets";
	progress_message(&oss,method_name);
#endif
	int error;
	int * node_list;
	float * dist_fact;

	my_node_set_ids = (int *) calloc(my_num_node_sets, sizeof(int));
	my_num_nodes_in_node_set = (int *) calloc(my_num_node_sets, sizeof(int));
	my_num_df_in_node_set = (int *) calloc(my_num_node_sets, sizeof(int));
	error = ex_get_node_set_ids (my_input_exoid, my_node_set_ids);
	for (int i=0; i<my_num_node_sets; i++)
	{
		error = ex_get_node_set_param (my_input_exoid, my_node_set_ids[i],
				&(my_num_nodes_in_node_set[i]), &(my_num_df_in_node_set[i]));
#ifdef DEBUG_OUTPUT
		oss << "Node set: " << my_node_set_ids[i];
		sub_progress_message(&oss);
		oss << "Number of nodes in set: " << my_num_nodes_in_node_set[i];
		sub_sub_progress_message(&oss);
		oss << "Number of distribution factors in set: " << my_num_df_in_node_set[i];
		sub_sub_progress_message(&oss);
#endif
		node_list = (int *) calloc(my_num_nodes_in_node_set[i], sizeof(int));
		dist_fact = (float *) calloc(my_num_nodes_in_node_set[i], sizeof(float));
		error = ex_get_node_set (my_input_exoid, my_node_set_ids[i], node_list);
		if (my_num_df_in_node_set[i] > 0)
		{
			error = ex_get_node_set_dist_fact (my_input_exoid, my_node_set_ids[i], dist_fact);
		}
        my_node_set_node_lists.insert(pair<int,int*>(my_node_set_ids[i],node_list));
        my_node_set_dist_factors.insert(pair<int,float*>(my_node_set_ids[i],dist_fact));
	}
}

void
Mesh_Manager::import_side_sets()
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::import_side_sets()";
	stringstream oss;
	oss << "Reading side sets";
	progress_message(&oss,method_name);
#endif
	int error;
	int * elem_list;
	int * node_list;
	float * dist_fact;
	int * side_list;
	int * node_ctr_list;

	my_side_set_ids = (int *) calloc(my_num_side_sets, sizeof(int));
	my_num_elem_in_side_set = (int *) calloc(my_num_side_sets, sizeof(int));
	my_num_df_in_side_set = (int *) calloc(my_num_side_sets, sizeof(int));
	error = ex_get_side_set_ids (my_input_exoid, my_side_set_ids);

	for (int i=0; i<my_num_side_sets; i++)
	{
		error = ex_get_side_set_param (my_input_exoid, my_side_set_ids[i], &(my_num_elem_in_side_set[i]),
				&(my_num_df_in_side_set[i]));
#ifdef DEBUG_OUTPUT
		oss << "Side set " << my_side_set_ids[i];
		sub_progress_message(&oss);
		oss << "Number of elements: " << my_num_elem_in_side_set[i];
		sub_sub_progress_message(&oss);
#endif
		/* Note: The # of elements is same as # of sides! */
		elem_list = (int *) calloc(my_num_elem_in_side_set[i], sizeof(int));
		side_list = (int *) calloc(my_num_elem_in_side_set[i], sizeof(int));
		node_ctr_list = (int *) calloc(my_num_elem_in_side_set[i], sizeof(int));
		node_list = (int *) calloc(my_num_elem_in_side_set[i]*21, sizeof(int));
		dist_fact = (float *) calloc(my_num_df_in_side_set[i], sizeof(float));
		error = ex_get_side_set (my_input_exoid, my_side_set_ids[i], elem_list, side_list);
		error = ex_get_side_set_node_list (my_input_exoid, my_side_set_ids[i], node_ctr_list,
				node_list);
		if (my_num_df_in_side_set[i] > 0)
		{
			error = ex_get_side_set_dist_fact (my_input_exoid, my_side_set_ids[i], dist_fact);
		}
        my_side_set_node_lists.insert(pair<int,int*>(my_side_set_ids[i],node_list));
        my_side_set_elem_lists.insert(pair<int,int*>(my_side_set_ids[i],elem_list));
        my_side_set_node_ctr_lists.insert(pair<int,int*>(my_side_set_ids[i],node_ctr_list));
        my_side_set_dist_factors.insert(pair<int,float*>(my_side_set_ids[i],dist_fact));
        my_side_set_side_lists.insert(pair<int,int*>(my_side_set_ids[i],side_list));
	}
}

void
Mesh_Manager::initialize_output(char * title)
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_output()";
	stringstream oss;
	oss << "Writing output file: " << my_output_file_name;
	progress_message(&oss,method_name);
#endif

	int error;

	my_CPU_word_size = 0;
	my_IO_word_size = 0;
	/* create EXODUS II file */
	my_output_exoid = ex_create (my_output_file_name,EX_CLOBBER,&my_CPU_word_size, &my_IO_word_size);

	error = ex_put_init (my_output_exoid, title, my_num_dim, my_num_nodes, my_num_elem, my_num_elem_blk, my_num_node_sets, my_num_side_sets);

    write_coordinates();
    write_elem_map();
    write_elem_blocks();
    write_elem_connectivities();
    write_node_sets();
    write_side_sets();
    write_qa_records();
    write_variable_names();

    // Truth table for element variables?
    //	/* write element variable truth table */
    //	truth_tab = (int *) calloc ((num_elem_blk*num_ele_vars), sizeof(int));
    //	k = 0;
    //	for (i=0; i<num_elem_blk; i++) {
    //		for (j=0; j<num_ele_vars; j++) {
    //			truth_tab[k++] = 1;
    //		}
    //	}
    //	error = ex_put_elem_var_tab (my_output_exoid, num_elem_blk, num_ele_vars, truth_tab); free (truth_tab);

    my_output_is_initialized = true;
}

void
Mesh_Manager::write_coordinates()
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_coordinates()";
	stringstream oss;
	oss << "Writing coordinates";
	progress_message(&oss,method_name);
#endif
	int error;
	error = ex_put_coord (my_output_exoid, my_x, my_y, my_z);
	char * coord_names[3];
	coord_names[0] = (char*)"xcoor";
	coord_names[1] = (char*)"ycoor";
	coord_names[2] = (char*)"zcoor";
	error = ex_put_coord_names (my_output_exoid, coord_names);
}

void
Mesh_Manager::write_elem_map()
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_elem_map()";
	stringstream oss;
	oss << "Writing element map";
	progress_message(&oss,method_name);
#endif
	int error;
	error = ex_put_map (my_output_exoid, my_elem_map);
}

void
Mesh_Manager::write_elem_blocks()
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_elem_blocks()";
	stringstream oss;
	oss << "Writing element blocks";
	progress_message(&oss,method_name);
#endif
	int error;
	for(int i=0; i<my_num_elem_blk;++i)
	{
		string elem_type_str = my_elem_types.find(my_block_ids[i])->second;
        char * elem_type = const_cast<char *>(elem_type_str.c_str());
		error = ex_put_elem_block (my_output_exoid, my_block_ids[i], elem_type, my_num_elem_in_block[i], my_num_nodes_per_elem[i], 0);
		/* write element block properties */
		//	prop_names[0] = (char*)"TOP"; prop_names[1] = (char*)"RIGHT"; error = ex_put_prop_names(my_output_exoid,EX_ELEM_BLOCK,2,prop_names);
		//	error = ex_put_prop(my_output_exoid, EX_ELEM_BLOCK, ebids[0], (char*)"TOP", 1);
		//	float attrib[0] = 6.14159;
		//	error = ex_put_elem_attr (my_output_exoid, my_block_ids[i], attrib);
		//	error = ex_put_elem_attr (my_output_exoid, my_block_ids[i], attrib);
	}
}

void
Mesh_Manager::write_elem_connectivities()
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_elem_connectivity()";
	stringstream oss;
	oss << "Writing element connectivity";
	progress_message(&oss,method_name);
#endif
	int error;
	int * connect;

	for(int i=0; i<my_num_elem_blk;++i)
	{
		connect = my_connectivities.find(my_block_ids[i])->second;
		error = ex_put_elem_conn (my_output_exoid, my_block_ids[i], connect);
	}
}

void
Mesh_Manager::write_node_sets()
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_node_sets()";
	stringstream oss;
	oss << "Writing node sets";
	progress_message(&oss,method_name);
#endif

	int error;
	for(int i=0; i<my_num_node_sets;++i)
	{
		error = ex_put_node_set_param (my_output_exoid, my_node_set_ids[i], my_num_nodes_in_node_set[i], my_num_df_in_node_set[i]);
		error = ex_put_node_set (my_output_exoid, my_node_set_ids[i], my_node_set_node_lists.find(my_node_set_ids[i])->second);
		error = ex_put_node_set_dist_fact (my_output_exoid, my_node_set_ids[i], my_node_set_dist_factors.find(my_node_set_ids[i])->second);
	    // error = ex_put_prop(my_output_exoid, EX_NODE_SET, 20, (char*)"FACE", 4);
		// prop_array[0] = 1000; prop_array[1] = 2000;
		// error = ex_put_prop_array(my_output_exoid, EX_NODE_SET, (char*)"VELOCITY", prop_array);
	}
}

void
Mesh_Manager::write_side_sets()
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_side_sets()";
	stringstream oss;
	oss << "Writing side sets";
	progress_message(&oss,method_name);
#endif
	int error;
	for(int i=0; i<my_num_side_sets;++i)
	{
		error = ex_put_side_set_param (my_output_exoid, my_side_set_ids[i], my_num_elem_in_side_set[i], my_num_df_in_side_set[i]);
		error = ex_put_side_set (my_output_exoid, my_side_set_ids[i], my_side_set_elem_lists.find(my_side_set_ids[i])->second, my_side_set_side_lists.find(my_side_set_ids[i])->second);
		error = ex_put_side_set_dist_fact (my_output_exoid, my_side_set_ids[i], my_side_set_dist_factors.find(my_side_set_ids[i])->second);
		//error = ex_put_prop(my_output_exoid, EX_SIDE_SET, 30, (char*)"COLOR", 100);
	}
}

void
Mesh_Manager::write_qa_records()
{
	stringstream oss;
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_qa_records()";
	oss << "Writing Q & A records";
	progress_message(&oss,method_name);
#endif
	int error;
	char *qa_record[1][4];
	int num_qa_rec = 1;

	qa_record[0][0] = (char*)"ANLYSIS CODE NAME GOES HERE";
	qa_record[0][1] = (char*)"code_name";

	struct tm *current;
	time_t now;
	time(&now);
	current = localtime(&now);
	int year = current->tm_year + 1900;
	int month = current->tm_mon + 1;
	oss << "  Date:    " << month << "/" <<  current->tm_mday << "/" << year << endl;
	char* date = const_cast<char *>(oss.str().c_str());
	oss.str(""); oss.clear();
	oss << current->tm_hour << ":" << current->tm_min << ":" << current->tm_sec;
	char* time = const_cast<char *>(oss.str().c_str());
	oss.str(""); oss.clear();

	qa_record[0][2] = date;
	qa_record[0][3] = time;

	error = ex_put_qa (my_output_exoid, num_qa_rec, qa_record);
}

void
Mesh_Manager::write_variable_names()
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_variable_names()";
	stringstream oss;
	oss << "Writing variable names";
	progress_message(&oss,method_name);
#endif
	int error;
	int num_global_variables = my_global_variable_names.size();
	char* global_var_names[num_global_variables];
	for(int i=0; i<num_global_variables;++i)
	{
        global_var_names[i] = (char *)my_global_variable_names[i];
	}
	error = ex_put_var_param (my_output_exoid, (char*)"g", num_global_variables);
	error = ex_put_var_names (my_output_exoid, (char*)"g", num_global_variables, global_var_names);
	int num_nodal_variables = my_nodal_variable_names.size();
	char* nodal_var_names[num_nodal_variables];
	for(int i=0; i<num_nodal_variables;++i)
	{
        nodal_var_names[i] = (char *)my_nodal_variable_names[i];
	}
	error = ex_put_var_param (my_output_exoid, (char*)"n", num_nodal_variables);
	error = ex_put_var_names (my_output_exoid, (char*)"n", num_nodal_variables, nodal_var_names);
	int num_element_variables = my_element_variable_names.size();
	char* ele_var_names[num_element_variables];
	for(int i=0; i<num_element_variables;++i)
	{
        ele_var_names[i] = (char *)my_element_variable_names[i];
	}
	error = ex_put_var_param (my_output_exoid, (char*)"e", num_element_variables);
	error = ex_put_var_names (my_output_exoid, (char*)"e", num_element_variables, ele_var_names);
}

void
Mesh_Manager::write_time_step_info(int time_step_num, float time_value)
{
	stringstream oss;
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_time_step_info()";
	oss << "Writing time step info to output";
	progress_message(&oss,method_name);
#endif
	int error;

	if(!my_output_is_initialized)
	{
		oss << "Output file is not initialized, can't write time step info to file: " << my_output_file_name;
		error_message(&oss);
	}
	error = ex_put_time (my_output_exoid, time_step_num, &time_value);
}

void
Mesh_Manager::update_output()
{
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::update_output()";
	stringstream oss;
	oss << "Updating output";
	progress_message(&oss,method_name);
#endif
	int error;
	/* update the data file; this should be done at the end of every time step * to ensure that no data is lost if the analysis dies */
	error = ex_update (my_output_exoid);
}

void
Mesh_Manager::write_global_variables_to_output(int time_step, float time_value, float * global_var_vals)
{
	stringstream oss;
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_global_variables_to_output()";
	oss << "Writing global variables to output";
	progress_message(&oss,method_name);
#endif
	int error;

	if(!my_output_is_initialized)
	{
		oss << "Output file is not initialized, can't write global variables to file: " << my_output_file_name;
		error_message(&oss);
	}
	int num_glo_vars = my_global_variable_names.size();
	error = ex_put_glob_vars (my_output_exoid, time_step, num_glo_vars, global_var_vals);
}

void
Mesh_Manager::write_nodal_variable_to_output(int time_step, float time_value, float * nodal_var_vals, int node_var_index)
{
	stringstream oss;
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_nodal_variable_to_output()";
	oss << "Writing nodal variable to output";
	progress_message(&oss,method_name);
#endif
	int error;

	if(!my_output_is_initialized)
	{
		oss << "Output file is not initialized, can't write nodal variables to file: " << my_output_file_name;
		error_message(&oss);
	}
	error = ex_put_nodal_var (my_output_exoid, time_step, node_var_index, my_num_nodes, nodal_var_vals);
}

void
Mesh_Manager::write_element_variable_to_output(int time_step, float time_value, float * elem_var_vals, int ele_var_index, int block_index)
{
	stringstream oss;
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_element_variable_to_output()";
	oss << "Writing element variable to output";
	progress_message(&oss,method_name);
#endif
	int error;
	if(!my_output_is_initialized)
	{
		oss << "Output file is not initialized, can't write variables to file: " << my_output_file_name;
		error_message(&oss);
	}
	error = ex_put_elem_var (my_output_exoid, time_step, ele_var_index, my_block_ids[block_index],my_num_elem_in_block[block_index], elem_var_vals);
}

void
Mesh_Manager::close_output_file()
{
	int error;
	error = ex_close (my_output_exoid);
}

void
Mesh_Manager::gregs_output()
{
	int exoid, num_dim, num_nodes, num_elem, num_elem_blk;
	int num_elem_in_block[10], num_nodes_per_elem[10];
	int num_node_sets, num_sides, num_side_sets, error;
	int i, j, k, m, *elem_map, *connect;
	int node_list[100],elem_list[100],side_list[100];
	int ebids[10], ids[10];
	int num_sides_per_set[10], num_nodes_per_set[10], num_elem_per_set[10];
	int num_df_per_set[10];
	int df_ind[10], node_ind[10], elem_ind[10], side_ind[10];
	int num_qa_rec, num_info;
	int num_glo_vars, num_nod_vars, num_ele_vars;
	int *truth_tab;
	int whole_time_step, num_time_steps;
	int ndims, nvars, ngatts, recdim;
	int CPU_word_size,IO_word_size;
	int prop_array[2];
	float *glob_var_vals, *nodal_var_vals, *elem_var_vals;
	float time_value;
	float x[100], y[100], z[100], *dummy;
	float attrib[1], dist_fact[100];
	char *coord_names[3], *qa_record[2][4], *info[3], *var_names[3];
	char tmpstr[80];
	char *prop_names[2];
	dummy = 0; /* assign this so the Cray compiler doesnÕt complain */ /* Specify compute and i/o word size */
	CPU_word_size = 0;/* float or double */ IO_word_size = 0;/* use system default (4 bytes) */
	/* create EXODUS II file */
	exoid = ex_create ("/Users/dzturne1/Documents/dzturne1/Research/cj/problems/unit_2d/gregs_output.e",/* filename path */ EX_CLOBBER,/* create mode */ &CPU_word_size,/* CPU float word size in bytes */&IO_word_size);/* I/O float word size in bytes */ /* ncopts = NC_VERBOSE; */
	/* initialize file with parameters */
	num_dim = 3; num_nodes = 26;
	num_elem = 5; num_elem_blk = 5; num_node_sets = 2; num_side_sets = 5;
	error = ex_put_init (exoid, "This is a test", num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets);
	/* write nodal coordinates values and names to database */
	/* Quad #1 */ x[0] = 0.0; y[0] = 0.0; z[0] = 0.0; x[1] = 1.0; y[1] = 0.0; z[1] = 0.0; x[2] = 1.0; y[2] = 1.0; z[2] = 0.0; x[3] = 0.0; y[3] = 1.0; z[3] = 0.0;
	/* Quad #2 */ x[4] = 1.0; y[4] = 0.0; z[4] = 0.0; x[5] = 2.0; y[5] = 0.0; z[5] = 0.0; x[6] = 2.0; y[6] = 1.0; z[6] = 0.0; x[7] = 1.0; y[7] = 1.0; z[7] = 0.0;
	/* Hex #1 */ x[8] = 0.0; y[8] = 0.0; z[8] = 0.0; x[9] = 10.0; y[9] = 0.0; z[9] = 0.0; x[10] = 10.0; y[10] = 0.0; z[10] =-10.0; x[11] = 1.0; y[11] = 0.0; z[11] =-10.0; x[12] = 1.0; y[12] = 10.0; z[12] = 0.0; x[13] = 10.0; y[13] = 10.0; z[13] = 0.0; x[14] = 10.0; y[14] = 10.0; z[14] =-10.0; x[15] = 1.0; y[15] = 10.0; z[15] =-10.0;
	/* Tetra #1 */ x[16] = 0.0; y[16] = 0.0; z[16] = 0.0; x[17] = 1.0; y[17] = 0.0; z[17] = 5.0; x[18] = 10.0; y[18] = 0.0; z[18] = 2.0; x[19] = 7.0; y[19] = 5.0; z[19] = 3.0;
	/* Wedge #1 */ x[20] = 3.0; y[20] = 0.0; z[20] = 6.0; x[21] = 6.0; y[21] = 0.0; z[21] = 0.0; x[22] = 0.0; y[22] = 0.0; z[22] = 0.0; x[23] = 3.0; y[23] = 2.0; z[23] = 6.0; x[24] = 6.0; y[24] = 2.0; z[24] = 2.0; x[25] = 0.0; y[25] = 2.0; z[25] = 0.0;
	error = ex_put_coord (exoid, x, y, z);
	coord_names[0] = (char *)"xcoor"; coord_names[1] = (char *)"ycoor"; coord_names[2] = (char *)"zcoor";
	error = ex_put_coord_names (exoid, coord_names); /* write element order map */
	elem_map = (int *) calloc(num_elem, sizeof(int));
	for (i=1; i<=num_elem; i++) {
	elem_map[i-1] = i;
	} error = ex_put_map (exoid, elem_map); free (elem_map);
	/* write element block parameters */
	num_elem_in_block[0] = 1;
	num_elem_in_block[1] = 1;
	num_elem_in_block[2] = 1;
	num_elem_in_block[3] = 1;
	num_elem_in_block[4] = 1;
	num_nodes_per_elem[0] = 4;
	/* elements in block #1 are 4-node quads */
	num_nodes_per_elem[1] = 4;
	/* elements in block #2 are 4-node quads */
	num_nodes_per_elem[2] = 8;
	/* elements in block #3 are 8-node hexes */
	num_nodes_per_elem[3] = 4;
	/* elements in block #3 are 4-node tetras */
	num_nodes_per_elem[4] = 6;
	/* elements in block #3 are 6-node wedges */
	ebids[0] = 10; ebids[1] = 11; ebids[2] = 12; ebids[3] = 13; ebids[4] = 14;
	error = ex_put_elem_block (exoid, ebids[0], (char *)"QUAD", num_elem_in_block[0], num_nodes_per_elem[0], 1);
	error = ex_put_elem_block (exoid, ebids[1], (char *)"QUAD", num_elem_in_block[1], num_nodes_per_elem[1], 1);
	error = ex_put_elem_block (exoid, ebids[2], (char *)"HEX", num_elem_in_block[2], num_nodes_per_elem[2], 1);
	error = ex_put_elem_block (exoid, ebids[3], (char *)"TETRA", num_elem_in_block[3], num_nodes_per_elem[3], 1);
	error = ex_put_elem_block (exoid, ebids[4], (char *)"WEDGE", num_elem_in_block[4], num_nodes_per_elem[4], 1);
	/* write element block properties */
	prop_names[0] = (char *)"TOP"; prop_names[1] = (char *)"RIGHT"; error = ex_put_prop_names(exoid,EX_ELEM_BLOCK,2,prop_names);
	error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[0], (char *)"TOP", 1);
	error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[1], (char *)"TOP", 1);
	error = ex_put_prop(exoid, EX_ELEM_BLOCK, ebids[2], (char *)"RIGHT", 1);
	error = ex_put_prop(exoid, EX_ELEM_BLOCK,  ebids[3], (char *)"RIGHT", 1);
	error = ex_put_prop(exoid, EX_ELEM_BLOCK,ebids[4], (char *)"RIGHT", 1);
	/* write element connectivity */

	connect = (int *) calloc(8, sizeof(int)); connect[0] = 1; connect[1] = 2; connect[2] = 3; connect[3] = 4;
	error = ex_put_elem_conn (exoid, ebids[0], connect);
	connect[0] = 5; connect[1] = 6; connect[2] = 7; connect[3] = 8;
	error = ex_put_elem_conn (exoid, ebids[1], connect);
	connect[0] = 9; connect[1] = 10; connect[2] = 11; connect[3] = 12; connect[4] = 13; connect[5] = 14; connect[6] = 15; connect[7] = 16;
	error = ex_put_elem_conn (exoid, ebids[2], connect);
	connect[0] = 17; connect[1] = 18; connect[2] = 19; connect[3] = 20;
	error = ex_put_elem_conn (exoid, ebids[3], connect);
	connect[0] = 21; connect[1] = 22; connect[2] = 23; connect[3] = 24; connect[4] = 25; connect[5] = 26;
	error = ex_put_elem_conn (exoid, ebids[4], connect);
	free (connect);
	/* write element block attributes */
	attrib[0] = 3.14159; error = ex_put_elem_attr (exoid, ebids[0], attrib);
	attrib[0] = 6.14159; error = ex_put_elem_attr (exoid, ebids[1], attrib);
	error = ex_put_elem_attr (exoid, ebids[2], attrib); error = ex_put_elem_attr (exoid, ebids[3], attrib); error = ex_put_elem_attr (exoid, ebids[4], attrib);
	/* write individual node sets */ error = ex_put_node_set_param (exoid, 20, 5, 5);
	node_list[0] = 100; node_list[1] = 101; node_list[2] = 102; node_list[3] = 103; node_list[4] = 104;
	dist_fact[0] = 1.0; dist_fact[1] = 2.0; dist_fact[2] = 3.0; dist_fact[3] = 4.0; dist_fact[4] = 5.0;
	error = ex_put_node_set (exoid, 20, node_list); error = ex_put_node_set_dist_fact (exoid, 20, dist_fact);
	error = ex_put_node_set_param (exoid, 21, 3, 3);
	node_list[0] = 200; node_list[1] = 201; node_list[2] = 202;
	dist_fact[0] = 1.1; dist_fact[1] = 2.1; dist_fact[2] = 3.1;
	error = ex_put_node_set (exoid, 21, node_list); error = ex_put_node_set_dist_fact (exoid, 21, dist_fact);
	error = ex_put_prop(exoid, EX_NODE_SET, 20, (char *)"FACE", 4); error = ex_put_prop(exoid, EX_NODE_SET, 21, (char *)"FACE", 5);
	prop_array[0] = 1000; prop_array[1] = 2000;
	error = ex_put_prop_array(exoid, EX_NODE_SET, (char *)"VELOCITY", prop_array);

	/* write individual side sets */
	/* side set #1 - quad */
	error = ex_put_side_set_param (exoid, 30, 2, 4);
	elem_list[0] = 2; elem_list[1] = 2;
	side_list[0] = 4; side_list[1] = 2;
	dist_fact[0] = 30.0; dist_fact[1] = 30.1; dist_fact[2] = 30.2; dist_fact[3] = 30.3;
	error = ex_put_side_set (exoid, 30, elem_list, side_list); error = ex_put_side_set_dist_fact (exoid, 30, dist_fact); /* side set #2 - quad, spanning 2 elements */ error = ex_put_side_set_param (exoid, 31, 2, 4); elem_list[0] = 1; elem_list[1] = 2;
	side_list[0] = 2; side_list[1] = 3;
	dist_fact[0] = 31.0; dist_fact[1] = 31.1; dist_fact[2] = 31.2; dist_fact[3] = 31.3;
	error = ex_put_side_set (exoid, 31, elem_list, side_list); error = ex_put_side_set_dist_fact (exoid, 31, dist_fact); /* side set #3 - hex */ error = ex_put_side_set_param (exoid, 32, 7, 0);
	elem_list[0] = 3; elem_list[2] = 3; elem_list[4] = 3; elem_list[6] = 3;
	side_list[0] = 5; side_list[2] = 3; side_list[4] = 4; side_list[6] = 6;
	elem_list[1] = 3; elem_list[3] = 3; elem_list[5] = 3;
	side_list[1] = 3; side_list[3] = 2; side_list[5] = 1;
	error = ex_put_side_set (exoid, 32,	elem_list, side_list);
	/* side set #4 - tetras */
	error = ex_put_side_set_param (exoid, 33, 4, 0);
	elem_list[0] = 4; elem_list[1] = 4; elem_list[2] = 4; elem_list[3] = 4;
	side_list[0] = 1; side_list[1] = 2; side_list[2] = 3; side_list[3] = 4;
	error = ex_put_side_set (exoid, 33, 	elem_list, side_list);
	/* side set #5 - wedges */
	error = ex_put_side_set_param (exoid, 34, 5, 0);

	elem_list[0] = 5; elem_list[2] = 5; elem_list[4] = 5;
	side_list[0] = 1; side_list[2] = 3; side_list[4] = 5;
	elem_list[1] = 5; elem_list[3] = 5;
	side_list[1] = 2; side_list[3] = 4;
	error = ex_put_side_set (exoid, 34,elem_list, side_list);


	error = ex_put_prop(exoid, EX_SIDE_SET, 30, (char *)"COLOR", 100);
	error = ex_put_prop(exoid, EX_SIDE_SET, 31, (char *)"COLOR", 101); /* write QA records */
	num_qa_rec = 2;
	qa_record[0][0] = (char *)"TESTWT";
	qa_record[0][1] = (char *)"testwt";
	qa_record[0][2] = (char *)"07/07/93";
	qa_record[0][3] = (char *)"15:41:33";
	qa_record[1][0]  = (char *)"FASTQ";
	qa_record[1][1] = (char *)"fastq";
	qa_record[1][2] = (char *)"07/07/93";
	qa_record[1][3] = (char *)"16:41:33";

	error = ex_put_qa (exoid, num_qa_rec, qa_record); /* write information records */
	num_info = 3;
	info[0] = (char *)"This is the first information record."; info[1] = (char *)"This is the second information record."; info[2] = (char *)"This is the third information record.";
	error = ex_put_info (exoid, num_info, info); /* write results variables parameters and names */
	num_glo_vars = 1;
	var_names[0] = (char *)"glo_vars";
	error = ex_put_var_param (exoid, "g", num_glo_vars);
	error = ex_put_var_names (exoid, "g",num_glo_vars, var_names);
	num_nod_vars = 2;
	var_names[0] = (char *)"nod_var0";
	var_names[1] = (char *)"nod_var1";
	error = ex_put_var_param (exoid, "n", num_nod_vars);
	error = ex_put_var_names (exoid, "n",num_nod_vars, var_names);
	num_ele_vars = 3;
	var_names[0] = (char *)"ele_var0";
	var_names[1] = (char *)"ele_var1";
	var_names[2] = (char *)"ele_var2";
	error = ex_put_var_param (exoid, "e", num_ele_vars);
	error = ex_put_var_names (exoid, "e",num_ele_vars, var_names);


	/* write element variable truth table */ truth_tab = (int *) calloc ((num_elem_blk*num_ele_vars), sizeof(int));
	k = 0; for (i=0; i<num_elem_blk; i++) {
	for (j=0; j<num_ele_vars; j++) {
	}
	truth_tab[k++] = 1;
	} error = ex_put_elem_var_tab (exoid, num_elem_blk, num_ele_vars, truth_tab); free (truth_tab);
	/* for each time step, write the analysis results; * the code below fills the arrays glob_var_vals,
	D-8
	* nodal_var_vals, and elem_var_vals with values for debugging purposes; * obviously the analysis code will populate these arrays */
	whole_time_step = 1; num_time_steps = 10;
	glob_var_vals = (float *) calloc (num_glo_vars, CPU_word_size);
	nodal_var_vals = (float *) calloc (num_nodes, CPU_word_size);
	elem_var_vals = (float *) calloc (4, CPU_word_size);
	for (i=0; i<num_time_steps; i++) {
		time_value = (float)(i+1)/100.; /* write time value */
		error = ex_put_time (exoid, whole_time_step, &time_value); /* write global variables */
		for (j=0; j<num_glo_vars; j++) {
			glob_var_vals[j] = (float)(j+2) * time_value;
		}
		error = ex_put_glob_vars (exoid, whole_time_step, num_glo_vars, glob_var_vals);
		/* write nodal variables */
		for (k=1; k<=num_nod_vars; k++) {
			for (j=0; j<num_nodes; j++) {
				nodal_var_vals[j] = (float)k + ((float)(j+1) * time_value);
			}
			error = ex_put_nodal_var (exoid, whole_time_step, k, num_nodes, nodal_var_vals);
		} /* write element variables */
		for (k=1; k<=num_ele_vars; k++) {
			for (j=0; j<num_elem_blk; j++) {
				for (m=0; m<num_elem_in_block[j]; m++) {
					elem_var_vals[m] = (float)(k+1) + (float)(j+2) + ((float)(m+1)*time_value);
				}
				error = ex_put_elem_var (exoid, whole_time_step, k, ebids[j],	num_elem_in_block[j], elem_var_vals);
			}
		}
		whole_time_step++;
		/* update the data file; this should be done at the end of every time step * to ensure that no data is lost if the analysis dies */
		error = ex_update (exoid);
	}
	free(glob_var_vals); free(nodal_var_vals); free(elem_var_vals);

	/* close the EXODUS files */
	error = ex_close (exoid);
}

