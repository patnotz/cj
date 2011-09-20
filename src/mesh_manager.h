#ifndef MESH_MANAGER_H_
#define MESH_MANAGER_H_

#include <vector>
#include <map>
#include "exodusII.h"
#include <stk_mesh.h>

class Mesh_Manager {
public:
	Mesh_Manager(const char * input_file_name, const char * output_file_name);
	~Mesh_Manager();

	void print_field_info(stk::mesh::STK_Mesh * const mesh);

	void populate_bogus_scalar_field(stk::mesh::STK_Mesh * const mesh,stk::mesh::ScalarFieldType & field, const stk::mesh::EntityRank & rank);
	void populate_bogus_vector_field(stk::mesh::STK_Mesh * const mesh,stk::mesh::VectorFieldType & field, const stk::mesh::EntityRank & rank);

	void populate_mesh_coordinates(stk::mesh::STK_Mesh * const mesh);
	void populate_mesh_elements(stk::mesh::STK_Mesh * const mesh);
	void initialize_mesh_parts_and_commit(stk::mesh::STK_Mesh * const mesh);
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
			double coord[]);
	bool verify_coordinates_field(const stk::mesh::STK_Mesh & mesh);

	void read_mesh();
	void initialize_output(const char * title, stk::mesh::STK_Mesh & mesh);
	void insert_global_var_name(const char * name) {
		my_global_variable_names.push_back(name);
	}
	void insert_nodal_var_name(const char * name)
	{
		for(int i=0;i<my_nodal_variable_names.size();++i)
		{
			if(my_nodal_variable_names[i]==name)
				return;
		}
		my_nodal_variable_names.push_back(name);
	}
    int get_nodal_variable_index(stk::mesh::STK_Mesh * const mesh,const std::string & name, const std::string & component);
    int get_element_variable_index(stk::mesh::STK_Mesh * const mesh,const std::string & name, const std::string & component);

    std::vector<std::string> get_field_names(stk::mesh::STK_Mesh * const mesh, stk::mesh::EntityRank entity_rank, unsigned field_rank);

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
	const int num_blocks() {
		return my_num_elem_blk;
	}
	const int num_elem_in_block(const int block_id) {
		if (block_id < my_num_elem_blk)
			return my_num_elem_in_block[block_id];
		else
			return 0;
	}

	void write_time_step_info(const int & time_step_num, const float & time_value);
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
	void write_nodal_scalar(const int & time_step,const float & time_value, stk::mesh::STK_Mesh & mesh,const stk::mesh::ScalarFieldType & field);
	void write_nodal_vector(const int & time_step,const float & time_value, stk::mesh::STK_Mesh & mesh,const stk::mesh::VectorFieldType & field);
	void write_element_scalar(const int & time_step,const float & time_value, stk::mesh::STK_Mesh & mesh,const stk::mesh::ScalarFieldType & field);
	void write_element_vector(const int & time_step,const float & time_value, stk::mesh::STK_Mesh & mesh,const stk::mesh::VectorFieldType & field);

	bool bucket_blocks_contain_block(const stk::mesh::Bucket & bucket, const std::string & name);

private:
	void initialize_read();
	void import_nodes();
	void import_elem_map();
	void import_blocks();
	void import_connectivities();
	void print_connectivity(const int & block_id);
	void import_node_sets();
	void import_side_sets();

	void write_coordinates(const stk::mesh::STK_Mesh & mesh);
	void write_nodal_vector(const stk::mesh::STK_Mesh & mesh,const stk::mesh::VectorFieldType & field);
	void write_elem_map();
	void write_elem_blocks(const stk::mesh::STK_Mesh & mesh);
	void write_elem_connectivities();
	void write_node_sets();
	void write_side_sets();
	void write_qa_records();
	void write_variable_names(stk::mesh::STK_Mesh * const  mesh);

	const char * my_input_file_name;
	const char * my_output_file_name;

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

	bool my_output_initialized;
	bool my_input_initialized;
};

#endif /* MESH_MANAGER_H_ */
