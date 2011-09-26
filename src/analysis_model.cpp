#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <sstream>
#include <time.h>

#include <analysis_model.h>
#include <physics.h>
#include <messages.h>
#include <log.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <Shards_CellTopologyData.h>
#include <Shards_BasicTopologies.hpp>

#include <exodusII.h>
#include <netcdf.h>


using namespace std;

Analysis_Model_Factory::Analysis_Model_Factory(){}

Teuchos::RCP<Analysis_Model> Analysis_Model_Factory::create(const char * input_file_name, const char * output_file_name, const std::string & name)
{
  using Teuchos::rcp;
  return rcp(new Analysis_Model(input_file_name, output_file_name, name));
}

Analysis_Model::Analysis_Model(const char * input_file_name, const char * output_file_name, const std::string & name) :
    my_input_file_name(input_file_name),
    my_output_file_name(output_file_name),
    my_input_initialized(false),
    my_output_initialized(false),
    my_name(name)
{
  stringstream oss1;
  oss1 << my_input_file_name;
  my_input_file_name_str = oss1.str();
  stringstream oss2;
  oss2 << my_output_file_name;
  my_output_file_name_str = oss2.str();
}

int Analysis_Model::get_nodal_variable_index(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh,
  const std::string & name, const std::string & component)
{
  string comp = "";
  if (component != "") comp = "_" + component;
  // this needs to be improved... right now the exodusII file does not differente between scalars and vectors for the variable names
  // so we have to concatinate the list like this
  // build up the list
  std::vector<std::string> nodal_scalar_fields = get_field_names(mesh,
    mesh->my_node_rank, 0);
  std::vector<std::string> nodal_vector_fields = get_field_names(mesh,
    mesh->my_node_rank, 1);
  int num_nodal_scalar_variables = nodal_scalar_fields.size();
  int num_nodal_vector_variables = nodal_vector_fields.size();
  int num_nodal_variables = num_nodal_scalar_variables
      + num_nodal_vector_variables;
  std::vector<std::string> nodal_var_names;
  for (int i = 0; i < num_nodal_scalar_variables; ++i)
  {
    nodal_var_names.push_back(nodal_scalar_fields[i]);
  }
  for (int i = 0; i < num_nodal_vector_variables; ++i)
  {
    nodal_var_names.push_back(nodal_vector_fields[i] + "_x");
    nodal_var_names.push_back(nodal_vector_fields[i] + "_y");
    if (mesh->my_spatial_dimension > 2) nodal_var_names.push_back(
      nodal_vector_fields[i] + "_z");
  }

  // rip through the list looking for the specified field which is broken down by component
  for (int i = 0; i < nodal_var_names.size(); ++i)
  {
    if (nodal_var_names[i] == name + comp) //FIXME this is case sensitive
    return i + 1;
  }
  return -1;
}

int Analysis_Model::get_element_variable_index(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh,
  const std::string & name, const std::string & component)
{
  string comp = "";
  if (component != "") comp = "_" + component;
  // this needs to be improved... right now the exodusII file does not differente between scalars and vectors for the variable names
  // so we have to concatinate the list like this
  // build up the list
  std::vector<std::string> element_scalar_fields = get_field_names(mesh,
    mesh->my_elem_rank, 0);
  std::vector<std::string> element_vector_fields = get_field_names(mesh,
    mesh->my_elem_rank, 1);
  int num_element_scalar_variables = element_scalar_fields.size();
  int num_element_vector_variables = element_vector_fields.size();
  int num_element_variables = num_element_scalar_variables
      + num_element_vector_variables;
  std::vector<std::string> element_var_names;
  for (int i = 0; i < num_element_scalar_variables; ++i)
  {
    element_var_names.push_back(element_scalar_fields[i]);
  }
  for (int i = 0; i < num_element_vector_variables; ++i)
  {
    element_var_names.push_back(element_vector_fields[i] + "_x");
    element_var_names.push_back(element_vector_fields[i] + "_y");
    if (mesh->my_spatial_dimension > 2) element_var_names.push_back(
      element_vector_fields[i] + "_z");
  }

  // rip through the list looking for the specified field which is broken down by component
  for (int i = 0; i < element_var_names.size(); ++i)
  {
    if (element_var_names[i] == name + comp) //FIXME this is case sensitive
    return i + 1;
  }
  return -1;
}

void Analysis_Model::map_node_coordinates(stk::mesh::EntityId node_id,
  double coord[]) const
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::map_node_coordinates()";
  oss << "Adding the coordinates field for node " << node_id
      << " to STK mesh field my_coordinates";
  progress_message(&oss, method_name);
#endif
  if (node_id < 1 || node_id > my_num_nodes)
  {
    oss << "map_node_coordinates(): ERROR, node (" << node_id
        << ") must be greater than 0 or less than " << my_num_nodes
        << std::endl;
    error_message(&oss);
    return;
  }
  const unsigned index = node_id - 1;
  coord[0] = my_x[index];
  coord[1] = my_y[index];
  if (my_num_dim > 2) coord[2] = my_z[index];
}

void Analysis_Model::map_node_ids(const int block, const int ele,
  stk::mesh::EntityId node_ids[], const string & elem_type,
  const int * connectivity) const
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::map_node_ids()";
  oss << "Converting the ExodusII connectivity for element " << ele
      << " of block " << block << " to STK mesh. Element type is " << elem_type;
  progress_message(&oss, method_name);
#endif
  if (ele < 0 || ele >= my_num_elem_in_block[block])
  {
    oss << "map_node_ids(): ERROR, element (" << ele
        << ") must be greater than 0 or less than "
        << my_num_elem_in_block[block] << std::endl;
    error_message(&oss);
    return;
  }

  if (elem_type == "QUAD4")
  {
    const unsigned base = (ele) * 4;
    node_ids[0] = connectivity[base + 3];
    node_ids[1] = connectivity[base + 0];
    node_ids[2] = connectivity[base + 1];
    node_ids[3] = connectivity[base + 2];
  }
  else if (elem_type == "HEX8")
  {
    const unsigned base = (ele) * 8;
    node_ids[0] = connectivity[base + 1];
    node_ids[1] = connectivity[base + 5];
    node_ids[2] = connectivity[base + 6];
    node_ids[3] = connectivity[base + 2];
    node_ids[4] = connectivity[base + 0];
    node_ids[5] = connectivity[base + 4];
    node_ids[6] = connectivity[base + 7];
    node_ids[7] = connectivity[base + 3];
  }
  else if (elem_type == "TETRA4" || elem_type == "TETRA")
  {
    const unsigned base = (ele) * 4;
    node_ids[0] = connectivity[base + 0];
    node_ids[1] = connectivity[base + 1];
    node_ids[2] = connectivity[base + 2];
    node_ids[3] = connectivity[base + 3];
  }
  else if (elem_type == "TRI3")
  {
    const unsigned base = (ele) * 3;
    node_ids[0] = connectivity[base + 0];
    node_ids[1] = connectivity[base + 1];
    node_ids[2] = connectivity[base + 2];
  }
  else if (elem_type == "SPHERE")
  {
    // For a sphere mesh or particle mesh the connectivity is just one node (itself)
    const unsigned base = (ele);
    node_ids[0] = connectivity[base + 0];
  }
  else
  {
    oss << "map_node_ids() does not recognize element type: " << elem_type;
    error_message(&oss);
    exit(1);
  }
}

void Analysis_Model::read_mesh()
{
  stringstream oss;
  int error;
  float version;

  int CPU_word_size = 0;
  int IO_word_size = 0;
  /*open exodus II file */

  const char * file_name = input_file_name().c_str();

  my_input_exoid =
      ex_open(file_name,EX_READ,&CPU_word_size,&IO_word_size,&version);
  if (my_input_exoid < 0)
  {
    oss << "Reading mesh failure: " << file_name;
    error_message(&oss);
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
  my_input_initialized = true;
}

void Analysis_Model::initialize_read()
{
  string method_name = "Analysis_Model::initialize_read()";
  int error;
  stringstream oss;

  const char * file_name = input_file_name().c_str();

#ifdef DEBUG_OUTPUT
  oss << "Reading file: " << file_name;
  progress_message(&oss, method_name);
#endif
  /* read database parameters */
  char title[MAX_LINE_LENGTH + 1];

error  = ex_get_init(my_input_exoid, title, &my_num_dim, &my_num_nodes, &my_num_elem,
    &my_num_elem_blk, &my_num_node_sets, &my_num_side_sets);

log()
    << endl
    << "  --------------------------------------------------------------------------------------"
      << endl
      << "    Analysis model: " << my_name
      << endl
      << "    Title: "
      << title
      << endl
      << "    Input file: " << my_input_file_name_str << "    Output file: " << my_output_file_name_str
      << endl
      << "    Spatial dimension: "
      << my_num_dim
      << endl
      << "    "
      << my_num_nodes
      << " node(s)   "
      << my_num_elem
      << " element(s)   "
      << my_num_elem_blk
      << " block(s)   "
      << my_num_node_sets
      << " node set(s)   "
      << my_num_side_sets
      << " side set(s)"
      << endl
      << "  --------------------------------------------------------------------------------------"
      << endl;
}

void Analysis_Model::import_nodes()
{
  int error;
  /* read nodal coordinates values and names from database */
  my_x = new float[my_num_nodes];
  my_y = new float[my_num_nodes];
  if (my_num_dim >= 3) my_z = new float[my_num_nodes];
  else
    my_z = 0;
  error = ex_get_coord(my_input_exoid, my_x, my_y, my_z);
}

void Analysis_Model::import_elem_map()
{
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::import_elem_map()";
  stringstream oss;
  oss << "Reading element map";
  progress_message(&oss, method_name);
#endif

  int error;
  /* read element order map */
  my_elem_map = new int[my_num_elem]; //FIXME: free me later
  error = ex_get_map(my_input_exoid, my_elem_map);
}

void Analysis_Model::import_blocks()
{
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::import_blocks()";
  stringstream oss;
  oss << "Reading blocks";
  progress_message(&oss, method_name);
#endif

  char elem_type[MAX_STR_LENGTH + 1];int
  error;

  /* read element block parameters */
  my_block_ids = new int[my_num_elem_blk]; // FIXME: free me later
  my_num_elem_in_block = new int[my_num_elem_blk]; // FIXME: free me later
  my_num_nodes_per_elem = new int[my_num_elem_blk]; // FIXME: free me later
  my_num_attr = new int[my_num_elem_blk]; // FIXME: free me later
  error = ex_get_elem_blk_ids(my_input_exoid, my_block_ids);
  for (int i = 0; i < my_num_elem_blk; i++)
  {
    error = ex_get_elem_block(my_input_exoid, my_block_ids[i], elem_type,
      &(my_num_elem_in_block[i]), &(my_num_nodes_per_elem[i]),
      &(my_num_attr[i]));

    my_elem_types.insert(pair<int, string>(my_block_ids[i], elem_type));

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

void Analysis_Model::import_connectivities()
{
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::import_connectivity()";
  stringstream oss;
#endif

  int error;
  int * connectivity;

  /* read element connectivity */
  for (int i = 0; i < my_num_elem_blk; i++)
  {
#ifdef DEBUG_OUTPUT
    oss << "Reading connectivity for block " << my_block_ids[i];
    progress_message(&oss, method_name);
#endif
    connectivity = new int[my_num_nodes_per_elem[i] * my_num_elem_in_block[i]]; // FIXME: free me later
error    = ex_get_elem_conn(my_input_exoid, my_block_ids[i], connectivity);
    my_connectivities.insert(pair<int, int*>(my_block_ids[i], connectivity));
  }
}

void Analysis_Model::print_connectivity(const int & block_id)
{
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::print_connectivity()";
  stringstream oss;
  oss << "Block ID: " << block_id;
  progress_message(&oss, method_name);
#endif
  if (my_connectivities.find(block_id) == my_connectivities.end()) return;

  int * connectivity = my_connectivities.find(block_id)->second;

  int num_elem = my_num_elem_in_block[block_id - 1];
  int num_nodes_per_elem = my_num_nodes_per_elem[block_id - 1];
  for (int i = 0; i < num_elem; ++i)
  {
    log() << "Element " << i << ": ";
    for (int j = 0; j < num_nodes_per_elem; ++j)
    {
      const int index = i * num_nodes_per_elem + j;
      log() << connectivity[index] << " ";
    }
    log() << endl;
  }
}

void Analysis_Model::import_node_sets()
{
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::import_node_sets()";
  stringstream oss;
  oss << "Reading node sets";
  progress_message(&oss, method_name);
#endif
  int error;
  int * node_list;
  float * dist_fact;

  my_node_set_ids = new int[my_num_node_sets]; // FIXME: free me later
  my_num_nodes_in_node_set = new int[my_num_node_sets]; // FIXME: free me later
  my_num_df_in_node_set = new int[my_num_node_sets]; // FIXME: free me later
  error = ex_get_node_set_ids(my_input_exoid, my_node_set_ids);
  for (int i = 0; i < my_num_node_sets; i++)
  {
    error = ex_get_node_set_param(my_input_exoid, my_node_set_ids[i],
      &(my_num_nodes_in_node_set[i]), &(my_num_df_in_node_set[i]));
#ifdef DEBUG_OUTPUT
    oss << "Node set: " << my_node_set_ids[i];
    sub_progress_message(&oss);
    oss << "Number of nodes in set: " << my_num_nodes_in_node_set[i];
    sub_sub_progress_message(&oss);
    oss << "Number of distribution factors in set: "
        << my_num_df_in_node_set[i];
    sub_sub_progress_message(&oss);
#endif
    node_list = new int[my_num_nodes_in_node_set[i]]; // FIXME: free me later
    dist_fact = new float[my_num_nodes_in_node_set[i]]; // FIXME: free me later
    error = ex_get_node_set(my_input_exoid, my_node_set_ids[i], node_list);
    if (my_num_df_in_node_set[i] > 0)
    {
      error = ex_get_node_set_dist_fact(my_input_exoid, my_node_set_ids[i],
        dist_fact);
    }
    my_node_set_node_lists.insert(
      pair<int, int*>(my_node_set_ids[i], node_list));
    my_node_set_dist_factors.insert(
      pair<int, float*>(my_node_set_ids[i], dist_fact));
  }
}

void Analysis_Model::import_side_sets()
{
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::import_side_sets()";
  stringstream oss;
  oss << "Reading side sets";
  progress_message(&oss, method_name);
#endif
  int error;
  int * elem_list;
  int * node_list;
  float * dist_fact;
  int * side_list;
  int * node_ctr_list;

  my_side_set_ids = new int[my_num_side_sets]; // FIXME: free me later
  my_num_elem_in_side_set = new int[my_num_side_sets]; // FIXME: free me later
  my_num_df_in_side_set = new int[my_num_side_sets]; // FIXME: free me later
  error = ex_get_side_set_ids(my_input_exoid, my_side_set_ids);

  for (int i = 0; i < my_num_side_sets; i++)
  {
    error = ex_get_side_set_param(my_input_exoid, my_side_set_ids[i],
      &(my_num_elem_in_side_set[i]), &(my_num_df_in_side_set[i]));
#ifdef DEBUG_OUTPUT
    oss << "Side set " << my_side_set_ids[i];
    sub_progress_message(&oss);
    oss << "Number of elements: " << my_num_elem_in_side_set[i];
    sub_sub_progress_message(&oss);
#endif
    /* Note: The # of elements is same as # of sides! */
    elem_list = new int[my_num_elem_in_side_set[i]]; // FIXME: free me later
    side_list = new int[my_num_elem_in_side_set[i]]; // FIXME: free me later
    node_ctr_list = new int[my_num_elem_in_side_set[i]]; // FIXME: free me later
    node_list = new int[my_num_elem_in_side_set[i] * 21]; // FIXME: free me later
dist_fact    = new float[my_num_df_in_side_set[i]]; // FIXME: free me later
    error = ex_get_side_set(my_input_exoid, my_side_set_ids[i], elem_list,
      side_list);
    error = ex_get_side_set_node_list(my_input_exoid, my_side_set_ids[i],
      node_ctr_list, node_list);
    if (my_num_df_in_side_set[i] > 0)
    {
      error = ex_get_side_set_dist_fact(my_input_exoid, my_side_set_ids[i],
        dist_fact);
    }
    my_side_set_node_lists.insert(
      pair<int, int*>(my_side_set_ids[i], node_list));
    my_side_set_elem_lists.insert(
      pair<int, int*>(my_side_set_ids[i], elem_list));
    my_side_set_node_ctr_lists.insert(
      pair<int, int*>(my_side_set_ids[i], node_ctr_list));
    my_side_set_dist_factors.insert(
      pair<int, float*>(my_side_set_ids[i], dist_fact));
    my_side_set_side_lists.insert(
      pair<int, int*>(my_side_set_ids[i], side_list));
  }
}

void Analysis_Model::initialize_output(const char * title,
  Teuchos::RCP<stk::mesh::STK_Mesh> const mesh)
{
  const char * file_name = output_file_name().c_str();

  #ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_output()";
  stringstream oss;

  oss << "Writing output file: " << file_name;
  progress_message(&oss, method_name);
#endif

  int error;

  int CPU_word_size = 0;
  int IO_word_size = 0;
  /* create EXODUS II file */
  my_output_exoid =
      ex_create (file_name,EX_CLOBBER,&CPU_word_size, &IO_word_size);

  error = ex_put_init(my_output_exoid, title, my_num_dim, my_num_nodes,
    my_num_elem, my_num_elem_blk, my_num_node_sets, my_num_side_sets);

  write_coordinates(mesh);
  //	const stk::mesh::VectorFieldType & coordinates_field = mesh.my_coordinates_field ;
  //    write_nodal_vector(mesh,coordinates_field);

  write_elem_map();
  write_elem_blocks();
  write_elem_connectivities();
  write_node_sets();
  write_side_sets();
  write_qa_records();
  write_variable_names(mesh);

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

  my_output_initialized = true;
}

void Analysis_Model::write_coordinates(Teuchos::RCP<stk::mesh::STK_Mesh> mesh)
{
  stringstream oss;

#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_coordinates()";
  oss << "Writing coordinates";
  progress_message(&oss, method_name);
#endif
  int error;

  float x[my_num_nodes];
  float y[my_num_nodes];
  float z[my_num_nodes];

  bool result = true;

  const stk::mesh::VectorFieldType & coordinates_field =
      mesh->my_coordinates_field;
  const stk::mesh::BulkData & bulkData = mesh->my_bulkData;
  const std::vector<stk::mesh::Bucket*> & node_buckets =
      mesh->my_bulkData.buckets(mesh->my_node_rank);

  for (std::vector<stk::mesh::Bucket*>::const_iterator node_bucket_it =
      node_buckets.begin(); node_bucket_it != node_buckets.end();
      ++node_bucket_it)
      {
    const stk::mesh::Bucket & bucket = **node_bucket_it;

    stk::mesh::BucketArray<stk::mesh::VectorFieldType> coordinates_array(
      mesh->my_coordinates_field, bucket);
    const int num_sets_of_coords = coordinates_array.dimension(1); //this is linked to the bucket size
    for (int i = 0; i < num_sets_of_coords; ++i)
    {
      const unsigned node_id = bucket[i].identifier();
      const unsigned index = node_id - 1;
      x[index] = coordinates_array(0, i);
      y[index] = coordinates_array(1, i);
      if (my_num_dim > 2) z[index] = coordinates_array(2, i);
    }
  }
  error = ex_put_coord(my_output_exoid, x, y, z);
  char * coord_names[3];
  coord_names[0] = (char*) "xcoor";
  coord_names[1] = (char*) "ycoor";
  coord_names[2] = (char*) "zcoor";
  error = ex_put_coord_names(my_output_exoid, coord_names);
}

void Analysis_Model::write_nodal_scalar(const int & time_step,
  const float & time_value, Teuchos::RCP<stk::mesh::STK_Mesh> mesh,
  const stk::mesh::ScalarFieldType & field)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_nodal_scalar()";
  oss << "Writing nodal scalar: " << field.name();
  progress_message(&oss, method_name);
#endif
  int error;
  float values[my_num_nodes];
  const std::vector<stk::mesh::Bucket*> & node_buckets =
      mesh->my_bulkData.buckets(mesh->my_node_rank);

  for (std::vector<stk::mesh::Bucket*>::const_iterator node_bucket_it =
      node_buckets.begin(); node_bucket_it != node_buckets.end();
      ++node_bucket_it)
      {
    const stk::mesh::Bucket & bucket = **node_bucket_it;
    stk::mesh::BucketArray<stk::mesh::ScalarFieldType> field_array(field,
      bucket);
    const int num_nodes_in_bucket = field_array.dimension(0); //this is linked to the bucket size
    for (int i = 0; i < num_nodes_in_bucket; ++i)
    {
      const unsigned node_id = bucket[i].identifier();
      const unsigned index = node_id - 1;
      values[index] = field_array(i);
    }
  }
  const int var_index = get_nodal_variable_index(mesh, field.name(), "");
  write_nodal_variable_to_output(time_step, time_value, values, var_index);
}

void Analysis_Model::write_nodal_vector(const int & time_step,
  const float & time_value, Teuchos::RCP<stk::mesh::STK_Mesh> mesh,
  const stk::mesh::VectorFieldType & field)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_nodal_vector()";
  oss << "Writing nodal vector: " << field.name();
  progress_message(&oss, method_name);
#endif
  int error;
  float values[my_num_nodes];
  const std::vector<stk::mesh::Bucket*> & node_buckets =
      mesh->my_bulkData.buckets(mesh->my_node_rank);

  string components[3];
  components[0] = "x";
  components[1] = "y";
  components[2] = "z";
  for (int comp = 0; comp < mesh->my_spatial_dimension; ++comp)
  {
    for (std::vector<stk::mesh::Bucket*>::const_iterator node_bucket_it =
        node_buckets.begin(); node_bucket_it != node_buckets.end();
        ++node_bucket_it)
        {
      const stk::mesh::Bucket & bucket = **node_bucket_it;
      stk::mesh::BucketArray<stk::mesh::VectorFieldType> field_array(field,
        bucket);
      const int num_nodes_in_bucket = field_array.dimension(1); //this is linked to the bucket size
      for (int i = 0; i < num_nodes_in_bucket; ++i)
      {
        const unsigned node_id = bucket[i].identifier();
        const unsigned index = node_id - 1;
        values[index] = field_array(comp, i);
      }
    }
    const int var_index = get_nodal_variable_index(mesh, field.name(),
      components[comp]);
    write_nodal_variable_to_output(time_step, time_value, values, var_index);
  }

}

bool bucket_blocks_contain_block_named(const stk::mesh::Bucket & bucket,
  const std::string & name)
{
  stk::mesh::PartVector elem_parts;
  bucket.supersets(elem_parts);
  for (int i = 0; i < elem_parts.size(); ++i)
  {
    if (elem_parts[i]->name() == name) return true;
  }
  return false;
}

void Analysis_Model::write_element_scalar(const int & time_step,
  const float & time_value, Teuchos::RCP<stk::mesh::STK_Mesh> mesh,
  const stk::mesh::ScalarFieldType & field)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_element_scalar()";
  oss << "Writing element scalar: " << field.name();
  progress_message(&oss, method_name);
#endif
  int error;

  for (int blk = 0; blk < num_blocks(); blk++)
  {
    int array_index = 0; //reset the index for each block
    const int num_elements = num_elem_in_block(blk);
    float values[num_elements];
    const std::vector<stk::mesh::Bucket*> & elem_buckets =
        mesh->my_bulkData.buckets(mesh->my_elem_rank); // get all of the element buckets

    // TODO this could be done more efficiently, but I don't know how yet.
    // As it is we rip over all the buckets for each block and only keep the
    // buckets that are part of this block. we need to figure out how to assemble
    // elements based on the block_id (mesh part)

    for (std::vector<stk::mesh::Bucket*>::const_iterator elem_bucket_it =
        elem_buckets.begin(); elem_bucket_it != elem_buckets.end();
        ++elem_bucket_it)
        {
      const stk::mesh::Bucket & bucket = **elem_bucket_it;
      bool bucket_has_block = bucket_blocks_contain_block_named(bucket,
        mesh->my_parts[blk]->name());
      if (!bucket_has_block) continue;
      // Now populate the array of element values for this block tha come from this bucket
      stk::mesh::BucketArray<stk::mesh::ScalarFieldType> field_array(field,
        bucket);
      const int num_elems_in_bucket = field_array.dimension(0); //this is linked to the bucket size
      for (int i = 0; i < num_elems_in_bucket; ++i)
      {
        const unsigned elem_id = bucket[i].identifier();
        //const unsigned index = elem_id - 1;
        values[array_index] = field_array(i);
        array_index++;
      }
    }
    const int var_index = get_element_variable_index(mesh, field.name(), "");
    write_element_variable_to_output(time_step, time_value, values, var_index,
      blk);
  }
}

void Analysis_Model::write_element_vector(const int & time_step,
  const float & time_value, Teuchos::RCP<stk::mesh::STK_Mesh> mesh,
  const stk::mesh::VectorFieldType & field)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_element_vector()";
  oss << "Writing element vector: " << field.name();
  progress_message(&oss, method_name);
#endif
  int error;

  string components[3];
  components[0] = "x";
  components[1] = "y";
  components[2] = "z";
  for (int comp = 0; comp < mesh->my_spatial_dimension; ++comp)
  {
    for (int blk = 0; blk < num_blocks(); blk++)
    {
      int array_index = 0;
      const int num_elements = num_elem_in_block(blk);
      float values[num_elements];
      const std::vector<stk::mesh::Bucket*> & elem_buckets =
          mesh->my_bulkData.buckets(mesh->my_elem_rank); // get all of the element buckets

      // TODO this could be done more efficiently, but I don't know how yet.
      // As it is we rip over all the buckets for each block and only keep the
      // buckets that are part of this block. we need to figure out how to assemble
      // elements based on the block_id (mesh part)

      for (std::vector<stk::mesh::Bucket*>::const_iterator elem_bucket_it =
          elem_buckets.begin(); elem_bucket_it != elem_buckets.end();
          ++elem_bucket_it)
          {
        const stk::mesh::Bucket & bucket = **elem_bucket_it;

        bool bucket_has_block = bucket_blocks_contain_block_named(bucket,
          mesh->my_parts[blk]->name());
        if (!bucket_has_block) continue;
        // Now populate the array of element values for this block tha come from this bucket
        stk::mesh::BucketArray<stk::mesh::VectorFieldType> field_array(field,
          bucket);
        const int num_elems_in_bucket = field_array.dimension(1); //this is linked to the bucket size
        for (int i = 0; i < num_elems_in_bucket; ++i)
        {
          values[array_index] = field_array(comp, i);
          array_index++;
        }
      }
      const int var_index = get_element_variable_index(mesh, field.name(),
        components[comp]);
      write_element_variable_to_output(time_step, time_value, values, var_index,
        blk);
    }
  }
}

void Analysis_Model::write_elem_map()
{
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_elem_map()";
  stringstream oss;
  oss << "Writing element map";
  progress_message(&oss, method_name);
#endif
  int error;
  error = ex_put_map(my_output_exoid, my_elem_map);
}

void Analysis_Model::write_elem_blocks()
{
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_elem_blocks()";
  stringstream oss;
  oss << "Writing element blocks";
  progress_message(&oss, method_name);
#endif
  int error;

  for (int i = 0; i < my_num_elem_blk; ++i)
  {
    string elem_type_str = my_elem_types.find(my_block_ids[i])->second;
    char * elem_type = const_cast<char *>(elem_type_str.c_str());error
    = ex_put_elem_block(my_output_exoid, my_block_ids[i], elem_type,
      my_num_elem_in_block[i], my_num_nodes_per_elem[i], 0);
    /* write element block properties */
    //	prop_names[0] = (char*)"TOP"; prop_names[1] = (char*)"RIGHT"; error = ex_put_prop_names(my_output_exoid,EX_ELEM_BLOCK,2,prop_names);
    //	error = ex_put_prop(my_output_exoid, EX_ELEM_BLOCK, ebids[0], (char*)"TOP", 1);
    //	float attrib[0] = 6.14159;
    //	error = ex_put_elem_attr (my_output_exoid, my_block_ids[i], attrib);
    //	error = ex_put_elem_attr (my_output_exoid, my_block_ids[i], attrib);
  }
}

void Analysis_Model::write_elem_connectivities()
{
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_elem_connectivity()";
  stringstream oss;
  oss << "Writing element connectivity";
  progress_message(&oss, method_name);
#endif
  int error;
  int * connect;

  for (int i = 0; i < my_num_elem_blk; ++i)
  {
    connect = my_connectivities.find(my_block_ids[i])->second;
    error = ex_put_elem_conn(my_output_exoid, my_block_ids[i], connect);
  }
}

void Analysis_Model::write_node_sets()
{
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_node_sets()";
  stringstream oss;
  oss << "Writing node sets";
  progress_message(&oss, method_name);
#endif

  int error;
  for (int i = 0; i < my_num_node_sets; ++i)
  {
    error = ex_put_node_set_param(my_output_exoid, my_node_set_ids[i],
      my_num_nodes_in_node_set[i], my_num_df_in_node_set[i]);
    error = ex_put_node_set(my_output_exoid, my_node_set_ids[i],
      my_node_set_node_lists.find(my_node_set_ids[i])->second);
    error = ex_put_node_set_dist_fact(my_output_exoid, my_node_set_ids[i],
      my_node_set_dist_factors.find(my_node_set_ids[i])->second);
    // error = ex_put_prop(my_output_exoid, EX_NODE_SET, 20, (char*)"FACE", 4);
    // prop_array[0] = 1000; prop_array[1] = 2000;
    // error = ex_put_prop_array(my_output_exoid, EX_NODE_SET, (char*)"VELOCITY", prop_array);
  }
}

void Analysis_Model::write_side_sets()
{
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_side_sets()";
  stringstream oss;
  oss << "Writing side sets";
  progress_message(&oss, method_name);
#endif
  int error;
  for (int i = 0; i < my_num_side_sets; ++i)
  {
    error = ex_put_side_set_param(my_output_exoid, my_side_set_ids[i],
      my_num_elem_in_side_set[i], my_num_df_in_side_set[i]);
    error = ex_put_side_set(my_output_exoid, my_side_set_ids[i],
      my_side_set_elem_lists.find(my_side_set_ids[i])->second,
      my_side_set_side_lists.find(my_side_set_ids[i])->second);
    error = ex_put_side_set_dist_fact(my_output_exoid, my_side_set_ids[i],
      my_side_set_dist_factors.find(my_side_set_ids[i])->second);
    //error = ex_put_prop(my_output_exoid, EX_SIDE_SET, 30, (char*)"COLOR", 100);
  }
}

void Analysis_Model::write_qa_records()
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_qa_records()";
  oss << "Writing Q & A records";
  progress_message(&oss, method_name);
#endif
  int error;
  char *qa_record[1][4];
  int num_qa_rec = 1;

  qa_record[0][0] = (char*) "ANLYSIS CODE NAME GOES HERE";
  qa_record[0][1] = (char*) "code_name";

  struct tm *current;
  time_t now;
  time(&now);
  current = localtime(&now);
  int year = current->tm_year + 1900;
  int month = current->tm_mon + 1;
  oss << "  Date:    " << month << "/" << current->tm_mday << "/" << year
      << endl;
  char* date = const_cast<char *>(oss.str().c_str());oss
  .str("");
  oss.clear();
  oss << current->tm_hour << ":" << current->tm_min << ":" << current->tm_sec;
  char* time = const_cast<char *>(oss.str().c_str());oss
  .str("");
  oss.clear();

  qa_record[0][2] = date;
  qa_record[0][3] = time;

  error = ex_put_qa(my_output_exoid, num_qa_rec, qa_record);
}

void Analysis_Model::write_variable_names(Teuchos::RCP<stk::mesh::STK_Mesh> mesh)
{
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_variable_names()";
  stringstream oss;
  oss << "Writing variable names";
  progress_message(&oss, method_name);
#endif
  int error;
  // NODAL fields
  std::vector<std::string> nodal_scalar_fields = get_field_names(mesh,
    mesh->my_node_rank, 0);
  std::vector<std::string> nodal_vector_fields = get_field_names(mesh,
    mesh->my_node_rank, 1);
  int num_nodal_scalar_variables = nodal_scalar_fields.size();
  int num_nodal_vector_variables = nodal_vector_fields.size();
  int num_nodal_variables = num_nodal_scalar_variables
      + num_nodal_vector_variables * mesh->my_spatial_dimension;
  std::vector<std::string> nodal_string_names;
  char* nodal_var_names[num_nodal_variables];
  for (int i = 0; i < num_nodal_scalar_variables; ++i)
  {
    nodal_string_names.push_back(nodal_scalar_fields[i]);
  }
  for (int i = 0; i < num_nodal_vector_variables; ++i)
  {
    string x_comp = nodal_vector_fields[i] + "_x";
    string y_comp = nodal_vector_fields[i] + "_y";
    string z_comp = nodal_vector_fields[i] + "_z";
    nodal_string_names.push_back(x_comp);
    nodal_string_names.push_back(y_comp);
    if (mesh->my_spatial_dimension > 2) nodal_string_names.push_back(z_comp);
  }
  for (int i = 0; i < num_nodal_variables; ++i)
  {
    nodal_var_names[i] = (char*) (nodal_string_names[i].c_str());
  }

  error = ex_put_var_param(my_output_exoid, (char*) "n", num_nodal_variables);
  error = ex_put_var_names(my_output_exoid, (char*) "n", num_nodal_variables,
    nodal_var_names);
  // ELEMENT fields
  std::vector<std::string> element_scalar_fields = get_field_names(mesh,
    mesh->my_elem_rank, 0);
  std::vector<std::string> element_vector_fields = get_field_names(mesh,
    mesh->my_elem_rank, 1);
  int num_element_scalar_variables = element_scalar_fields.size();
  int num_element_vector_variables = element_vector_fields.size();
  int num_element_variables = num_element_scalar_variables
      + num_element_vector_variables * mesh->my_spatial_dimension;
  char* ele_var_names[num_element_variables];
  std::vector<std::string> element_string_names;

  for (int i = 0; i < num_element_scalar_variables; ++i)
  {
    element_string_names.push_back(element_scalar_fields[i]);

  }
  for (int i = 0; i < num_element_vector_variables; ++i)
  {
    string x_comp = element_vector_fields[i] + "_x";
    string y_comp = element_vector_fields[i] + "_y";
    string z_comp = element_vector_fields[i] + "_z";
    element_string_names.push_back(x_comp);
    element_string_names.push_back(y_comp);
    if (mesh->my_spatial_dimension > 2) element_string_names.push_back(z_comp);

  }
  for (int i = 0; i < num_element_variables; ++i)
  {
    ele_var_names[i] = (char*) (element_string_names[i].c_str());
  }
  error = ex_put_var_param(my_output_exoid, (char*) "e", num_element_variables);
  error = ex_put_var_names(my_output_exoid, (char*) "e", num_element_variables,
    ele_var_names);
}

void Analysis_Model::write_time_step_info(const int & time_step_num,
  const float & time_value)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_time_step_info()";
  oss << "Writing time step info to output";
  progress_message(&oss, method_name);
#endif
  int error;

  const char * file_name = output_file_name().c_str();

  if (!my_output_initialized)
  {
    oss
        << "Output file is not initialized, can't write time step info to file: "
        << file_name;
    error_message(&oss);
  }
  error = ex_put_time(my_output_exoid, time_step_num, &time_value);
}

void Analysis_Model::update_output()
{
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::update_output()";
  stringstream oss;
  oss << "Updating output";
  progress_message(&oss, method_name);
#endif
  int error;
  /* update the data file; this should be done at the end of every time step * to ensure that no data is lost if the analysis dies */
  error = ex_update(my_output_exoid);
}

void Analysis_Model::write_global_variables_to_output(const int & time_step,
  const float & time_value, const float * global_var_vals)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_global_variables_to_output()";
  oss << "Writing global variables to output";
  progress_message(&oss, method_name);
#endif
  int error;

  const char * file_name = output_file_name().c_str();

  if (!my_output_initialized)
  {
    oss
        << "Output file is not initialized, can't write global variables to file: "
        << file_name;
    error_message(&oss);
  }
  int num_glo_vars = my_global_variable_names.size();
  error = ex_put_glob_vars(my_output_exoid, time_step, num_glo_vars,
    global_var_vals);
}

void Analysis_Model::write_nodal_variable_to_output(const int & time_step,
  const float & time_value, const float * nodal_var_vals,
  const int & node_var_index)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_nodal_variable_to_output()";
  oss << "Writing nodal variable to output";
  progress_message(&oss, method_name);
#endif
  int error;

  const char * file_name = output_file_name().c_str();

  if (!my_output_initialized)
  {
    oss
        << "Output file is not initialized, can't write nodal variables to file: "
        << file_name;
    error_message(&oss);
  }
  error = ex_put_nodal_var(my_output_exoid, time_step, node_var_index,
    my_num_nodes, nodal_var_vals);
}

void Analysis_Model::write_element_variable_to_output(const int & time_step,
  const float & time_value, const float * elem_var_vals,
  const int & ele_var_index, const int & block_index)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Analysis_Model::write_element_variable_to_output()";
  oss << "Writing element variable to output";
  progress_message(&oss, method_name);
#endif
  int error;
  const char * file_name = output_file_name().c_str();
  if (!my_output_initialized)
  {
    oss << "Output file is not initialized, can't write variables to file: "
        << file_name;
    error_message(&oss);
  }
  error = ex_put_elem_var(my_output_exoid, time_step, ele_var_index,
    my_block_ids[block_index], my_num_elem_in_block[block_index],
    elem_var_vals);
}

void Analysis_Model::close_output_file()
{
  int error;
  error = ex_close(my_output_exoid);
}

template<class field_type>
bool gather_field_data(unsigned expected_num_rel, const field_type & field,
  const stk::mesh::Entity & entity,
  typename stk::mesh::FieldTraits<field_type>::data_type * dst,
  stk::mesh::EntityRank entity_rank, const int dim)
{
  typedef typename stk::mesh::FieldTraits<field_type>::data_type T;

  stk::mesh::PairIterRelation rel = entity.relations(entity_rank);

  bool result = expected_num_rel == (unsigned) rel.size();

  if (result)
  {
    // Iterate over field data for each related entity and copy data
    // into src for one entity at a time
    T * const dst_end = dst + dim * expected_num_rel;
    for (; dst < dst_end; ++rel, dst += dim)
    {
      const T* src = field_data(field, *rel->entity());
      if (!src)
      {
        break;
      }
      // FIXME:
      if (dim == 1) stk::Copy<1>(dst, src);
      if (dim == 2) stk::Copy<2>(dst, src);
      if (dim == 3) stk::Copy<3>(dst, src);
    }
    result = dst == dst_end;
  }
  return result;
}
