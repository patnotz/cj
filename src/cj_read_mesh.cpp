#include <iostream>
#include <../include/cj_read_mesh.h>

#include "exodusII.h"
#include "netcdf.h"

using namespace std;

void
read_mesh(const char* mesh_file_name)
{
  int error;
  float version;
  int IO_word_size = 0, CPU_word_size=0, exoid=0;

  cout << "Reading mesh file: " << mesh_file_name << endl;

  exoid = ex_open(mesh_file_name,EX_READ,&CPU_word_size,&IO_word_size,&version);
  if(exoid<0)
  {
	  cout << "Reading mesh failure." << endl;
	  exit(1);
  }
  cout << "Mesh read successfully" << endl;

  /* read database parameters */
  char title[MAX_LINE_LENGTH+1];
  int num_nodes,num_dim,num_elem,num_elem_blk, num_node_sets, num_side_sets;

  error = ex_get_init (exoid, title, &num_dim, &num_nodes, &num_elem,
  &num_elem_blk, &num_node_sets, &num_side_sets);

  cout << "-------------------------------------------------------------" << endl <<
		  "|                      Mesh Data                             " << endl <<
		  "|                                                            " << endl <<
		  "|  Title: " << title << endl <<
		  "|  " << endl <<
		  "|  Spatial dimension: " << num_dim << endl <<
		  "|  Number of nodes: " << num_nodes << endl <<
		  "|  Number of elements: " << num_elem << endl <<
		  "|  Number of blocks: " << num_elem_blk << endl <<
		  "|  Number of node sets: " << num_node_sets << endl <<
		  "|  Number of side sets: " << num_side_sets << endl << "| " << endl <<
		  "-------------------------------------------------------------" << endl;

//  float *x, *y, *z;
//
//  /* read nodal coordinates values and names from database */
//  x = (float *) calloc(num_nodes, sizeof(float));
//  y = (float *) calloc(num_nodes, sizeof(float));
//  if (num_dim >= 3)
//  z = (float *) calloc(num_nodes, sizeof(float));
//  else
//  D-11
//  z = 0;
//  error = ex_get_coord (exoid, x, y, z);
//  free (x);
//  free (y);
//  if (num_dim >= 3)
//  free (z);
//  for (i=0; i<num_dim; i++)
//  {
//  coord_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
//  }
//  error = ex_get_coord_names (exoid, coord_names);
//  for (i=0; i<num_dim; i++)
//  free(coord_names[i]);


}
