#include <iostream>
#include <string>
#include <sstream>
#include <../include/mesh_reader.h>
#include <../include/messages.h>


#include "exodusII.h"
#include "netcdf.h"

using namespace std;

void
Mesh_Reader::read_mesh()
{
  string method_name = "Mesh_Reader::read_mesh";
  int error;
  float version;
  int IO_word_size = 0, CPU_word_size=0, exoid=0;

  stringstream oss;
  oss << "Reading file: " << my_file_name;
  progress_message(&oss, method_name);

//  cout << "%  Reading mesh file: " << mesh_file_name << endl;

  /*open exodus II file */
  exoid = ex_open(my_file_name,EX_READ,&CPU_word_size,&IO_word_size,&version);
  if(exoid<0)
  {
	  oss << "Reading mesh failure.";
	  sub_progress_message(&oss);
	  exit(1);
  }
  oss << "Mesh read successfully";
  sub_progress_message(&oss);

  /* read database parameters */
  char title[MAX_LINE_LENGTH+1];
  int num_nodes,num_dim,num_elem,num_elem_blk, num_node_sets, num_side_sets;

  error = ex_get_init (exoid, title, &num_dim, &num_nodes, &num_elem,
  &num_elem_blk, &num_node_sets, &num_side_sets);

  cout << endl <<
		  "  ----------------------------------------------------------------------------" << endl <<
		  "    Title: " << title << endl <<
		  "    Spatial dimension: " << num_dim << endl <<
		  "    " << num_nodes << " node(s)   " << num_elem << " element(s)   " << num_elem_blk << " block(s)   " << num_node_sets << " node set(s)   " <<  num_side_sets << " side set(s)" << endl <<
		  "  ----------------------------------------------------------------------------" << endl;

  float *x, *y, *z;

  /* read nodal coordinates values and names from database */
  x = (float *) calloc(num_nodes, sizeof(float));
  y = (float *) calloc(num_nodes, sizeof(float));
  if (num_dim >= 3)
	  z = (float *) calloc(num_nodes, sizeof(float));
  else
	  z = 0;
  error = ex_get_coord (exoid, x, y, z);
//  for(int i=0; i < num_nodes; i++)
//  {
//	  cout << "Node " << i << " coords (" << x[i] <<"," << y[i] << ")" << endl;
//  }
  free (x);
  free (y);
  if (num_dim >= 3)
  free (z);

  int * ids;
  int * num_elem_in_block;
  int * num_attr;
  int * num_nodes_per_elem;
  char elem_type[MAX_STR_LENGTH+1];

  oss << "Reading blocks";
  progress_message(&oss,method_name);

  /* read element block parameters */
  ids = (int *) calloc(num_elem_blk, sizeof(int));
  num_elem_in_block = (int *) calloc(num_elem_blk, sizeof(int));
  num_nodes_per_elem = (int *) calloc(num_elem_blk, sizeof(int));
  num_attr = (int *) calloc(num_elem_blk, sizeof(int));
  error = ex_get_elem_blk_ids (exoid, ids);
  for (int i=0; i<num_elem_blk; i++)
  {
	  error = ex_get_elem_block (exoid, ids[i], elem_type,
			  &(num_elem_in_block[i]),
			  &(num_nodes_per_elem[i]), &(num_attr[i]));

	  oss << "Element block: " << i;
	  sub_progress_message(&oss);
	  oss << "Element type: " << elem_type;
	  sub_sub_progress_message(&oss);
	  oss << "Number of elements in block: " << num_elem_in_block[i];
	  sub_sub_progress_message(&oss);
	  oss << "Number of nodes per element: " << num_nodes_per_elem[i];
	  sub_sub_progress_message(&oss);
	  oss << "Number of attributes: " << num_attr[i];
	  sub_sub_progress_message(&oss); oss.clear();
  };

  /* read element block properties */

//  int num_props;
//  float fdum;
//  char * cdum;
//  char * prop_names[3];
//  int prop_value;

//  cout << "%  Reading blocks properties" << endl;
//
//  error = ex_inquire (exoid, EX_INQ_EB_PROP, &num_props, &fdum, cdum);
//  for (int i=0; i<num_props; i++)
//  {
//	  prop_names[i] = (char *) calloc((MAX_STR_LENGTH+1), sizeof(char));
//  }
//  error = ex_get_prop_names(exoid,EX_ELEM_BLOCK,prop_names);
//  for (int i=0; i<num_props; i++)
//  {
//	  cout << "   - Property name: " << prop_names[i] << endl;
//	  for (int j=0; j<num_elem_blk; j++)
//	  {
//		  error = ex_get_prop(exoid, EX_ELEM_BLOCK, ids[j], prop_names[i],
//				  &prop_value);
//		  cout << "      Block: " << j << " value: " << prop_value << endl;
//	  }
//  }
//  for (int i=0; i<num_props; i++)
//	  free(prop_names[i]);

  int * connect;
  /* read element connectivity */
  for (int i=0; i<num_elem_blk; i++)
  {
	  oss << "Reading connectivity for block " << i;
	  progress_message(&oss, method_name);

	  connect = (int *) calloc((num_nodes_per_elem[i] * num_elem_in_block[i]),
			  sizeof(int));
	  error = ex_get_elem_conn (exoid, ids[i], connect);
//	  for(int j=0;j<num_nodes_per_elem[i] * num_elem_in_block[i];++j)
//	  {
//		  cout << connect[j] << endl;
//	  }
	  free (connect);
  }

  free (ids);
  free (num_nodes_per_elem);
  free (num_attr);

  int num_nodes_in_set;
  int num_df_in_set;
  int * node_list;
  float * dist_fact;

  oss << "Reading node sets";
  progress_message(&oss,method_name);

  /* read individual node sets */
  ids = (int *) calloc(num_node_sets, sizeof(int));
  error = ex_get_node_set_ids (exoid, ids);
  for (int i=0; i<num_node_sets; i++)
  {
	  oss << "Node set: " << i;
	  sub_progress_message(&oss);
	  error = ex_get_node_set_param (exoid, ids[i],
			  &num_nodes_in_set, &num_df_in_set);
	  oss << "Number of nodes in set: " << num_nodes_in_set;
	  sub_sub_progress_message(&oss);

	  node_list = (int *) calloc(num_nodes_in_set, sizeof(int));
	  dist_fact = (float *) calloc(num_nodes_in_set, sizeof(float));
	  error = ex_get_node_set (exoid, ids[i], node_list);
	  if (num_df_in_set > 0)
	  {
		  error = ex_get_node_set_dist_fact (exoid, ids[i], dist_fact);
	  }
	  free (node_list);
	  free (dist_fact);
  }
  free(ids);

  int num_sides_in_set;
  int num_elem_in_set;
  int * elem_list;
  int * side_list;
  int * node_ctr_list;

  oss << "Reading side sets";
  progress_message(&oss,method_name);

  /* read individual side sets */
  ids = (int *) calloc(num_side_sets, sizeof(int));
  error = ex_get_side_set_ids (exoid, ids);
  for (int i=0; i<num_side_sets; i++)
  {
	  error = ex_get_side_set_param (exoid, ids[i], &num_sides_in_set,
			  &num_df_in_set);
	  oss << "Side set " << i;
	  sub_progress_message(&oss);
	  oss << "Number of elements: " << num_sides_in_set;
	  sub_sub_progress_message(&oss);
	  /* Note: The # of elements is same as # of sides! */
	  num_elem_in_set = num_sides_in_set;
	  elem_list = (int *) calloc(num_elem_in_set, sizeof(int));
	  side_list = (int *) calloc(num_sides_in_set, sizeof(int));
	  node_ctr_list = (int *) calloc(num_elem_in_set, sizeof(int));
	  node_list = (int *) calloc(num_elem_in_set*21, sizeof(int));
	  dist_fact = (float *) calloc(num_df_in_set, sizeof(float));
	  error = ex_get_side_set (exoid, ids[i], elem_list, side_list);
	  error = ex_get_side_set_node_list (exoid, ids[i], node_ctr_list,
			  node_list);
	  if (num_df_in_set > 0)
	  {
		  error = ex_get_side_set_dist_fact (exoid, ids[i], dist_fact);
	  }
	  free (elem_list);
	  free (side_list);
	  free (node_ctr_list);
	  free (node_list);
	  free (dist_fact);
  }



}
