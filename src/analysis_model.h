#ifndef ANALYSIS_MODEL_H_
#define ANALYSIS_MODEL_H_

#include <vector>
#include <map>
#include "exodusII.h"
#include <stk_mesh.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <iostream>

class Analysis_Model
{
public:
  Analysis_Model(const char * input_file_name, const char * output_file_name, const std::string & name);
  ~Analysis_Model(){ };

  void map_node_ids(const int block, const int ele,
    stk::mesh::EntityId node_ids[], const std::string & elem_type,
    const int * connectivity) const;
  void map_node_coordinates(stk::mesh::EntityId node_id, double coord[]) const;
  void read_mesh();
  void initialize_output(const char * title, Teuchos::RCP<stk::mesh::STK_Mesh> const mesh);

  void insert_global_var_name(const char * name){my_global_variable_names.push_back(name);}
  void insert_nodal_var_name(const char * name)
  {
    for (int i = 0; i < my_nodal_variable_names.size(); ++i)
    {
      if (my_nodal_variable_names[i] == name) return;
    }
    my_nodal_variable_names.push_back(name);
  }
  int get_nodal_variable_index(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh,
    const std::string & name, const std::string & component);
  int get_element_variable_index(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh,
    const std::string & name, const std::string & component);

  const bool input_initialized()const {return my_input_initialized;}
  const bool output_initialized()const {return my_output_initialized;}

  void insert_element_var_name(const char * name)
  {
    my_element_variable_names.push_back(name);
  }
  const int num_global_variables()
  {
    return my_global_variable_names.size();
  }
  const int num_nodal_variables()
  {
    return my_nodal_variable_names.size();
  }
  const int num_element_variables()
  {
    return my_element_variable_names.size();
  }
  void close_output_file();
  const int num_nodes()
  {
    return my_num_nodes;
  }
  const int num_elem()
  {
    return my_num_elem;
  }
  const int num_blocks() const
  {
    return my_num_elem_blk;
  }
  const int num_elem_in_block(const int block_id)
  const{
    if (block_id < my_num_elem_blk) return my_num_elem_in_block[block_id];
    else
      return 0;
  }

  std::string name(){return my_name;}

  void write_time_step_info(const int & time_step_num,
    const float & time_value);
  void update_output();
  void write_global_variables_to_output(const int & time_step,
    const float & time_value, const float * global_var_vals);

  void write_nodal_variable_to_output(const int & time_step,
    const float & time_value, const float * nodal_var_vals,
    const int & node_var_index);

  void write_element_variable_to_output(const int & time_step,
    const float & time_value, const float * elem_var_vals,
    const int & ele_var_index, const int & block_index);

  const int dimension()
  {
    return my_num_dim;
  }
  void write_nodal_scalar(const int & time_step, const float & time_value,
    Teuchos::RCP<stk::mesh::STK_Mesh> mesh, const stk::mesh::ScalarFieldType & field);
  void write_nodal_vector(const int & time_step, const float & time_value,
    Teuchos::RCP<stk::mesh::STK_Mesh> mesh, const stk::mesh::VectorFieldType & field);
  void write_element_scalar(const int & time_step, const float & time_value,
    Teuchos::RCP<stk::mesh::STK_Mesh> mesh, const stk::mesh::ScalarFieldType & field);
  void write_element_vector(const int & time_step, const float & time_value,
    Teuchos::RCP<stk::mesh::STK_Mesh> mesh, const stk::mesh::VectorFieldType & field);

  bool bucket_blocks_contain_block(const stk::mesh::Bucket & bucket,
    const std::string & name);

  const std::string input_file_name(){return my_input_file_name_str;}
  const std::string output_file_name(){return my_output_file_name_str;}

  const std::map<int, std::string> elem_types()const{return my_elem_types;}
  const int * block_ids()const{return my_block_ids;}
  const std::map<int, int*> connectivities()const{return my_connectivities;}
  const int * num_nodes_per_elem()const {return my_num_nodes_per_elem;}
  const int * elem_map()const {return my_elem_map;}

  const int num_dim() const{return my_num_dim;}
  const float * x_coord()const{return my_x;}
  const float * y_coord()const{return my_y;}
  const float * z_coord()const{return my_z;}

private:

  void initialize_read();
  void import_nodes();
  void import_elem_map();
  void import_blocks();
  void import_connectivities();
  void print_connectivity(const int & block_id);
  void import_node_sets();
  void import_side_sets();

  void write_coordinates(Teuchos::RCP<stk::mesh::STK_Mesh> mesh);
  void write_nodal_vector(const stk::mesh::STK_Mesh & mesh,
    const stk::mesh::VectorFieldType & field);
  void write_elem_map();
  void write_elem_blocks();
  void write_elem_connectivities();
  void write_node_sets();
  void write_side_sets();
  void write_qa_records();
  void write_variable_names(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh);

  const char * my_input_file_name;
  const char * my_output_file_name;
  std::string my_input_file_name_str;
  std::string my_output_file_name_str;

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

  std::string my_name;
};

class Analysis_Model_Factory
{
public:
  Analysis_Model_Factory();
  virtual ~Analysis_Model_Factory(){}
  virtual Teuchos::RCP<Analysis_Model> create(const char * input_file_name, const char * output_file_name, const std::string & name);

private:
  Analysis_Model_Factory(const Analysis_Model_Factory&);
  Analysis_Model_Factory& operator=(const Analysis_Model_Factory&);

protected:

};

#endif /* ANALYSIS_MODEL_H_ */
