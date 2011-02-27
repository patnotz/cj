/*
 * mesh_manager.h
 *
 *  Created on: Jan 30, 2011
 *      Author: dzturne1
 */
#ifndef MESH_MANAGER_H_
#define MESH_MANAGER_H_

#include <vector>
#include <map>

#include "exodusII.h"

#include <stk_mesh.h>
#include <log.h>

class Mesh_Manager {
public:
	Mesh_Manager(const char * input_file_name, const char * output_file_name);
	~Mesh_Manager();

	void populate_STK_mesh(stk::mesh::STK_Mesh * const mesh);
	stk::mesh::Part * const
	part_pointer(stk::mesh::STK_Mesh * const mesh, const std::string & elem_type, const std::string & name);
	void map_node_ids(
			const int block,
			const int ele,
			stk::mesh::EntityId node_ids[],
			const std::string & elem_type,
			const int * connectivity);
	void map_node_coordinates(
			stk::mesh::EntityId node_id,
			double coord[],
			const int spatial_dim);
	bool verify_coordinates_field(Log & log, const stk::mesh::STK_Mesh & mesh);

	void read_mesh(Log & log);
	void initialize_output(const char * title);
	void insert_global_var_name(const char * name) {
		my_global_variable_names.push_back(name);
	}
	void insert_nodal_var_name(const char * name) {
		my_nodal_variable_names.push_back(name);
	}
	void insert_element_var_name(const char * name) {
		my_element_variable_names.push_back(name);
	}
	const int num_global_variables() {
		return my_global_variable_names.size();
	}
	const int num_nodal_variables() {
		return my_nodal_variable_names.size();
	}
	const int num_element_variables() {
		return my_element_variable_names.size();
	}
	void close_output_file();
	const int num_nodes() {
		return my_num_nodes;
	}
	const int num_elem() {
		return my_num_elem;
	}
	const int cpu_word_size() {
		return my_CPU_word_size;
	}
	const int io_word_size() {
		return my_IO_word_size;
	}
	const int num_blocks() {
		return my_num_elem_blk;
	}
	const int num_elem_in_block(const int block_id) {
		if (block_id < my_num_elem_blk)
			return my_num_elem_in_block[block_id];
		else
			return 0;
	}

	void
			write_time_step_info(const int & time_step_num, const float & time_value);
	void update_output();
	void write_global_variables_to_output(
			const int & time_step,
			const float & time_value,
			const float * global_var_vals);
	void write_nodal_variable_to_output(
			const int & time_step,
			const float & time_value,
			const float * nodal_var_vals,
			const int & node_var_index);
	void write_element_variable_to_output(
			const int & time_step,
			const float & time_value,
			const float * elem_var_vals,
			const int & ele_var_index,
			const int & block_index);

	const int dimension() {
		return my_num_dim;
	}

private:
	void initialize_read(Log & log);
	void import_nodes();
	void import_elem_map();
	void import_blocks();
	void import_connectivities();
	void print_connectivity(Log & log, const int & block_id);
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

	const char * my_input_file_name;
	const char * my_output_file_name;

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
	std::map<int, std::string> my_elem_types;
	std::map<int, int*> my_connectivities;

	int * my_node_set_ids;
	int * my_num_nodes_in_node_set;
	int * my_num_df_in_node_set;
	std::map<int, int*> my_node_set_node_lists;
	std::map<int, float*> my_node_set_dist_factors;

	int * my_side_set_ids;
	int * my_num_elem_in_side_set;
	int * my_num_df_in_side_set;
	std::map<int, int*> my_side_set_node_lists;
	std::map<int, int*> my_side_set_elem_lists;
	std::map<int, int*> my_side_set_side_lists;
	std::map<int, int*> my_side_set_node_ctr_lists;
	std::map<int, float*> my_side_set_dist_factors;

	int my_output_exoid;

	std::vector<const char *> my_global_variable_names;
	std::vector<const char *> my_nodal_variable_names;
	std::vector<const char *> my_element_variable_names;
	//int my_num_edge_variables;
	//int my_num_face_variables;

	bool my_output_initialized;
	bool my_input_initialized;
};

#endif /* MESH_MANAGER_H_ */
