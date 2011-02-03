/*
 * mesh_manager.h
 *
 *  Created on: Jan 30, 2011
 *      Author: dzturne1
 */
#include <iostream>
#include <vector>
#include <map>

#include "exodusII.h"

#ifndef MESH_MANAGER_H_
#define MESH_MANAGER_H_

using namespace std;

class Mesh_Manager
{
//private:
//    void insert_var_name(string & name, std::vector<string> * vector);
public:
	void set_input_file(char * file_name){my_input_file_name = file_name;}
	void set_output_file(char * file_name){my_output_file_name = file_name;}
    void read_mesh();
    void initialize_output(char * title);
    void insert_global_var_name(char *  name){my_global_variable_names.push_back(name);}
    void insert_nodal_var_name(char * name){my_nodal_variable_names.push_back(name);}
    void insert_element_var_name(char * name){my_element_variable_names.push_back(name);}
    char * get_nodal_variable_name(int index){if(index<my_nodal_variable_names.size())return my_nodal_variable_names[index]; else return (char *)"";};
    char * get_element_variable_name(int index){if(index<my_element_variable_names.size())return my_element_variable_names[index]; else return (char *)"";};
    const int num_global_variables(){return my_global_variable_names.size();}
    const int num_nodal_variables(){return my_nodal_variable_names.size();}
    const int num_element_variables(){return my_element_variable_names.size();}
    void close_output_file();
    const int num_nodes(){return my_num_nodes;}
    const int num_elem(){return my_num_elem;}
    const int cpu_word_size(){return my_CPU_word_size;}
    const int io_word_size(){return my_IO_word_size;}
    const int num_blocks(){return my_num_elem_blk;}
    const int num_elem_in_block(int block_id){if(block_id<my_num_elem_blk)return my_num_elem_in_block[block_id];else return 0;}

    void write_time_step_info(int time_step_num, float time_value);
    void update_output();
    void write_global_variables_to_output(int time_step, float time_value, float * global_var_vals);
    void write_nodal_variable_to_output(int time_step, float time_value, float * nodal_var_vals, int node_var_index);
    void write_element_variable_to_output(int time_step, float time_value, float * elem_var_vals, int ele_var_index,int block_index);

    void gregs_output();

private:
    void initialize_read();
    void import_nodes();
    void import_elem_map();
    void import_blocks();
    void import_connectivities();
    void print_connectivity(int block_id);
    void import_node_sets();
    void import_side_sets();

    void write_coordinates();
    void write_elem_map();
    void write_elem_blocks();
    void write_elem_connectivities();
    void write_node_sets();
    void write_side_sets();
    void write_qa_records();
    void write_variable_names();


    char * my_input_file_name;
	char * my_output_file_name;

	int my_IO_word_size;
	int my_CPU_word_size;

	int my_input_exoid;
	int my_num_nodes;
	int my_num_dim;
	int my_num_elem;
	int my_num_elem_blk;
	int my_num_node_sets;
	int my_num_side_sets;

	float *my_x, *my_y, *my_z;

	int * my_elem_map;

	int * my_block_ids;
	int * my_num_elem_in_block;
	int * my_num_attr;
	int * my_num_nodes_per_elem;
	map<int,string> my_elem_types;
    map<int,int*> my_connectivities;

    int * my_node_set_ids;
    int * my_num_nodes_in_node_set;
    int * my_num_df_in_node_set;
    map<int,int*> my_node_set_node_lists;
    map<int,float*> my_node_set_dist_factors;

    int * my_side_set_ids;
    int * my_num_elem_in_side_set;
    int * my_num_df_in_side_set;
    map<int,int*> my_side_set_node_lists;
    map<int,int*> my_side_set_elem_lists;
    map<int,int*> my_side_set_side_lists;
    map<int,int*> my_side_set_node_ctr_lists;
    map<int,float*> my_side_set_dist_factors;

    int my_output_exoid;

    std::vector<char *> my_global_variable_names;
    std::vector<char *> my_nodal_variable_names;
    std::vector<char *> my_element_variable_names;
    //int my_num_edge_variables;
    //int my_num_face_variables;

    bool my_output_is_initialized;
};

#endif /* MESH_MANAGER_H_ */
