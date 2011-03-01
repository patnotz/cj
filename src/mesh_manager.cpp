/*
 * mesh_manager.cpp
 *
 *  Created on: Jan 30, 2011
 *      Author: dzturne1
 */

#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <time.h>

#include <mesh_manager.h>
#include <messages.h>
#include <log.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/FEMInterface.hpp>

#include <Shards_CellTopologyData.h>

#include <exodusII.h>
#include <netcdf.h>

using stk::mesh::fem::NODE_RANK;

using namespace std;

Mesh_Manager::Mesh_Manager(const char * input_file_name, const char * output_file_name)
:  my_input_file_name(input_file_name),
   my_output_file_name(output_file_name),
   my_input_initialized(false),
   my_output_initialized(false)
{}

Mesh_Manager::~Mesh_Manager()
{}

void
Mesh_Manager::populate_STK_mesh(stk::mesh::STK_Mesh * const mesh)
{
	stringstream oss;
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::populate_STK_mesh()";
	oss << "Populating STK mesh";
	progress_message(oss, method_name);
#endif
	if(!my_input_initialized)
	{
		oss << "Output not initialized yet.";
		error_message(&oss);
		exit(1);
	}
#ifdef DEBUG_OUTPUT
	oss << "Populating the stk PartVector";
	sub_progress_message(&oss);
#endif
	// populate the list of block parts
	for(int i = 0 ; i < my_num_elem_blk; ++i)
	{
		const string elem_type = my_elem_types.find(my_block_ids[i])->second;
		oss << "block_" << my_block_ids[i] << "_" << elem_type;
		const string part_name = oss.str(); oss.str(""); oss.clear();

		stk::mesh::Part * part_ptr = part_pointer(mesh, elem_type, part_name);
		mesh->my_parts.push_back(part_ptr);
	}
#ifdef DEBUG_OUTPUT
	for(int i = 0 ; i < my_num_elem_blk; ++i)
	{
		oss << mesh->my_parts[i]->name();
		sub_sub_progress_message(&oss);
	}
#endif
    mesh->my_metaData.commit(); // this has to be called before we can modify the mesh
#ifdef DEBUG_OUTPUT
	oss << "Mesh committed, adding the elements";
	sub_progress_message(&oss);
#endif
	mesh->my_bulkData.modification_begin(); // Begin modifying the mesh
	int ele_map_index = 0;
	for(int block = 0 ; block < my_num_elem_blk; ++block)
	{
		int * connectivity = my_connectivities.find(my_block_ids[block])->second;
		const string elem_type = my_elem_types.find(my_block_ids[block])->second;
		for(int ele = 0 ; ele < my_num_elem_in_block[block]; ++ele)
		{
			stk::mesh::EntityId node_ids[ my_num_nodes_per_elem[block] ];
			stk::mesh::EntityId elem_id = my_elem_map[ele_map_index];

			// Note declare_element expects a cell topology
			// to have been attached to the part
			map_node_ids(block,ele,node_ids,elem_type,connectivity);
#ifdef DEBUG_OUTPUT
			oss << "ExodusII: ";
 			int base = 0;
			if(elem_type == "QUAD4")
			{
				base = ( ele ) * 4 ;
			}
			else if(elem_type == "HEX8")
			{
				base = ( ele ) * 8 ;
			}
			else if (elem_type == "TETRA")
			{
				base = ( ele ) * 4 ;
			}
			else if (elem_type == "TRI3")
			{
				base = ( ele ) * 3 ;
			}
			for(int node = 0; node< my_num_nodes_per_elem[block]; ++node)
			{
				oss << connectivity[base + node] << ",";
			}
			oss << " -> STK mesh: ";
			for(int node = 0; node< my_num_nodes_per_elem[block]; ++node)
			{
				oss << node_ids[node] << ",";
			}
			sub_sub_progress_message(&oss);
#endif
#ifdef DEBUG_OUTPUT
	oss << "Adding element id: " << elem_id << " in block " << mesh->my_parts[block]->name();
	sub_sub_progress_message(&oss);
#endif
		    stk::mesh::declare_element(mesh->my_bulkData,*mesh->my_parts[block],elem_id,node_ids);
			ele_map_index ++;
		}
	}
	// Done modifying the mesh.
	// Modifications on the local parallel process are communicated
	// among processes, verified for consistency, and changes to
	// parallel shared/ghosted mesh entities are synchronized.
	mesh->my_bulkData.modification_end();

#ifdef DEBUG_OUTPUT
	oss << "Populating the coordinates field";
	sub_progress_message(&oss);
#endif
    const std::vector<stk::mesh::Bucket*> & node_buckets =
    		mesh->my_bulkData.buckets(NODE_RANK);

    for ( std::vector<stk::mesh::Bucket*>::const_iterator
    		node_bucket_it = node_buckets.begin() ;
    		node_bucket_it != node_buckets.end() ; ++node_bucket_it )
    {
    	const stk::mesh::Bucket & bucket = **node_bucket_it;
#ifdef DEBUG_OUTPUT
	oss << "Coordinates field for bucket " << bucket.key();
	sub_sub_progress_message(&oss);
#endif
    	// Fill the nodal coordinates.
    	// Create a multidimensional array view of the
    	// nodal coordinates field data for this bucket of nodes.
    	// The array is two dimensional ( Cartesian X NumberNodes )
    	// and indexed by ( 0..2 , 0..NumerNodes-1 )
    	stk::mesh::BucketArray<stk::mesh::VectorFieldType> coordinates_array( mesh->my_coordinates_field, bucket );
    	const int num_sets_of_coords = coordinates_array.dimension(1);  //this is linked to the bucket size
    	for ( int i=0 ; i < num_sets_of_coords ; ++i )
    	{
    		const unsigned node_id = bucket[i].identifier();
    		map_node_coordinates(node_id,& coordinates_array(0,i));
    	}
    }
    //Now that x, y, and z are in the stk mesh field we can delete the temp data storage
    delete my_x;
    delete my_y;
    if(my_num_dim > 2)
    	delete my_z;
}

void
Mesh_Manager::map_node_coordinates( stk::mesh::EntityId node_id , double coord[])
{
	stringstream oss;
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::map_node_coordinates()";
	oss << "Adding the coordinates field for node " << node_id << " to STK mesh field my_coordinates";
	progress_message(&oss, method_name);
#endif
	if ( node_id < 1 || node_id > my_num_nodes) {
		oss << "map_node_coordinates(): ERROR, node ("
				<< node_id << ") must be greater than 0 or less than "<< my_num_nodes << std::endl;
		error_message(&oss);
		return;
	}
	const unsigned index = node_id - 1;
	coord[0] = my_x[index];
	coord[1] = my_y[index];
	if(my_num_dim > 2)
		coord[2] = my_z[index];
}

void
Mesh_Manager::map_node_ids(const int block, const int ele, stk::mesh::EntityId node_ids[], const string & elem_type, const int * connectivity)
{
	stringstream oss;
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::map_node_ids()";
	oss << "Converting the ExodusII connectivity for element " << ele << " of block " << block << " to STK mesh. Element type is " << elem_type;
	progress_message(&oss, method_name);
#endif
	if ( ele < 0 || ele >= my_num_elem_in_block[block]) {
		oss << "map_node_ids(): ERROR, element ("
				<< ele << ") must be greater than 0 or less than "<< my_num_elem_in_block[block] << std::endl;
		error_message(&oss);
		return;
	}

	if(elem_type == "QUAD4")
	{
		const unsigned base = ( ele ) * 4 ;
		node_ids[0] = connectivity[base + 3] ;
		node_ids[1] = connectivity[base + 0] ;
		node_ids[2] = connectivity[base + 1] ;
		node_ids[3] = connectivity[base + 2] ;
	}
	else if(elem_type == "HEX8")
	{
		const unsigned base = ( ele ) * 8 ;
		node_ids[0] = connectivity[base + 1] ;
		node_ids[1] = connectivity[base + 5] ;
		node_ids[2] = connectivity[base + 6] ;
		node_ids[3] = connectivity[base + 2] ;
		node_ids[4] = connectivity[base + 0] ;
		node_ids[5] = connectivity[base + 4] ;
		node_ids[6] = connectivity[base + 7] ;
		node_ids[7] = connectivity[base + 3] ;
	}
	else if (elem_type == "TETRA")
	{
		const unsigned base = ( ele ) * 4 ;
		node_ids[0] = connectivity[base + 0] ;
		node_ids[1] = connectivity[base + 1] ;
		node_ids[2] = connectivity[base + 2] ;
		node_ids[3] = connectivity[base + 3] ;
	}
	else if (elem_type == "TRI3")
	{
		const unsigned base = ( ele ) * 3 ;
		node_ids[0] = connectivity[base + 0] ;
		node_ids[1] = connectivity[base + 1] ;
		node_ids[2] = connectivity[base + 2] ;
	}
	else
	{
		oss << "map_node_ids() does not recognize element type: " << elem_type;
		error_message(&oss);
		exit(1);
	}
}

stk::mesh::Part * const
Mesh_Manager::part_pointer(stk::mesh::STK_Mesh * const mesh, const string & elem_type, const string & name)
{
	stringstream oss;
	stk::mesh::Part * part_ptr;

	if(elem_type == "QUAD4")
	{
		stk::mesh::Part & part = stk::mesh::declare_part<shards::Quadrilateral<4> >(mesh->my_metaData, name);
		part_ptr = &part;
	}
	else if(elem_type == "HEX8")
	{
		stk::mesh::Part & part = stk::mesh::declare_part<shards::Hexahedron<8> >(mesh->my_metaData, name);
		part_ptr = &part;
	}
	else if (elem_type == "TETRA")
	{
		stk::mesh::Part & part = stk::mesh::declare_part<shards::Tetrahedron<4> >(mesh->my_metaData, name);
		part_ptr = &part;
	}
	else if (elem_type == "TRI3")
	{
		stk::mesh::Part & part = stk::mesh::declare_part<shards::Triangle<3> >(mesh->my_metaData, name);
		part_ptr = &part;
	}
	else
	{
		oss << "part_pointer() does not recognize element type: " << elem_type;
		error_message(&oss);
		exit(1);
	}
	return part_ptr;
}

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

	  log() << endl <<
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
	my_x = new float[my_num_nodes];
	my_y = new float[my_num_nodes];
	if (my_num_dim >= 3)
		my_z = new float[my_num_nodes];
	else
		my_z = 0;
	error = ex_get_coord (my_input_exoid, my_x, my_y, my_z);
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
	my_elem_map = new int[my_num_elem]; //FIXME: free me later
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
	my_block_ids = new int[my_num_elem_blk]; // FIXME: free me later
	my_num_elem_in_block = new int[my_num_elem_blk];// FIXME: free me later
	my_num_nodes_per_elem = new int[my_num_elem_blk];// FIXME: free me later
	my_num_attr = new int[my_num_elem_blk];// FIXME: free me later
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
		connectivity = new int[my_num_nodes_per_elem[i] * my_num_elem_in_block[i]]; // FIXME: free me later
		error = ex_get_elem_conn (my_input_exoid, my_block_ids[i], connectivity);
        my_connectivities.insert(pair<int,int*>(my_block_ids[i],connectivity));
	}
}

void
Mesh_Manager::print_connectivity(const int & block_id)
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
		log() << "Element " << i << ": ";
		for(int j = 0;j < num_nodes_per_elem;++j)
		{
			const int index = i*num_nodes_per_elem + j;
			log() << connectivity[index] << " ";
		}
		log() << endl;
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

	my_node_set_ids = new int [my_num_node_sets];// FIXME: free me later
	my_num_nodes_in_node_set = new int [my_num_node_sets];// FIXME: free me later
	my_num_df_in_node_set = new int [my_num_node_sets]; // FIXME: free me later
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
		node_list = new int[my_num_nodes_in_node_set[i]]; // FIXME: free me later
		dist_fact = new float[my_num_nodes_in_node_set[i]];// FIXME: free me later
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

	my_side_set_ids = new int[my_num_side_sets]; // FIXME: free me later
	my_num_elem_in_side_set = new int [my_num_side_sets]; // FIXME: free me later
	my_num_df_in_side_set = new int [my_num_side_sets]; // FIXME: free me later
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
		elem_list = new int[my_num_elem_in_side_set[i]]; // FIXME: free me later
		side_list = new int[my_num_elem_in_side_set[i]];// FIXME: free me later
		node_ctr_list = new int [my_num_elem_in_side_set[i]]; // FIXME: free me later
		node_list = new int[my_num_elem_in_side_set[i]*21]; // FIXME: free me later
		dist_fact = new float [my_num_df_in_side_set[i]]; // FIXME: free me later
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
Mesh_Manager::initialize_output(const char * title, const stk::mesh::STK_Mesh & mesh)
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

    write_coordinates(mesh);
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

    my_output_initialized = true;
}

void
Mesh_Manager::write_coordinates(const stk::mesh::STK_Mesh & mesh)
{
	stringstream oss;

#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_coordinates()";
	oss << "Writing coordinates";
	progress_message(&oss,method_name);
#endif
	int error;
	//	error = ex_put_coord (my_output_exoid, my_x, my_y, my_z);
	//	char * coord_names[3];
	//	coord_names[0] = (char*)"xcoor";
	//	coord_names[1] = (char*)"ycoor";
	//	coord_names[2] = (char*)"zcoor";
	//	error = ex_put_coord_names (my_output_exoid, coord_names);

   float x[my_num_nodes];
   float y[my_num_nodes];
   float z[my_num_nodes];

	bool result = true;

	const stk::mesh::VectorFieldType & coordinates_field = mesh.my_coordinates_field ;
	const stk::mesh::BulkData & bulkData = mesh.my_bulkData ;

    const std::vector<stk::mesh::Bucket*> & node_buckets =
    		mesh.my_bulkData.buckets(NODE_RANK);

    for ( std::vector<stk::mesh::Bucket*>::const_iterator
    		node_bucket_it = node_buckets.begin() ;
    		node_bucket_it != node_buckets.end() ; ++node_bucket_it )
    {
    	const stk::mesh::Bucket & bucket = **node_bucket_it;

    	stk::mesh::BucketArray<stk::mesh::VectorFieldType> coordinates_array( mesh.my_coordinates_field, bucket );
    	const int num_sets_of_coords = coordinates_array.dimension(1);  //this is linked to the bucket size
    	for ( int i=0 ; i < num_sets_of_coords ; ++i )
    	{
    		const unsigned node_id = bucket[i].identifier();
    		const unsigned index = node_id - 1;
    		x[index] = coordinates_array(0,i);
    		y[index] = coordinates_array(1,i);
    		if(my_num_dim > 2)
    			z[index] = coordinates_array(2,i);
    	}
    }

    //	int error;
    error = ex_put_coord (my_output_exoid, x, y, z);
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
Mesh_Manager::write_time_step_info(const int & time_step_num, const float & time_value)
{
	stringstream oss;
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_time_step_info()";
	oss << "Writing time step info to output";
	progress_message(&oss,method_name);
#endif
	int error;

	if(!my_output_initialized)
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
Mesh_Manager::write_global_variables_to_output(const int & time_step, const float & time_value, const float * global_var_vals)
{
	stringstream oss;
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_global_variables_to_output()";
	oss << "Writing global variables to output";
	progress_message(&oss,method_name);
#endif
	int error;

	if(!my_output_initialized)
	{
		oss << "Output file is not initialized, can't write global variables to file: " << my_output_file_name;
		error_message(&oss);
	}
	int num_glo_vars = my_global_variable_names.size();
	error = ex_put_glob_vars (my_output_exoid, time_step, num_glo_vars, global_var_vals);
}

void
Mesh_Manager::write_nodal_variable_to_output(const int & time_step, const float & time_value, const float * nodal_var_vals, const int & node_var_index)
{
	stringstream oss;
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_nodal_variable_to_output()";
	oss << "Writing nodal variable to output";
	progress_message(&oss,method_name);
#endif
	int error;

	if(!my_output_initialized)
	{
		oss << "Output file is not initialized, can't write nodal variables to file: " << my_output_file_name;
		error_message(&oss);
	}
	error = ex_put_nodal_var (my_output_exoid, time_step, node_var_index, my_num_nodes, nodal_var_vals);
}

void
Mesh_Manager::write_element_variable_to_output(const int & time_step, const float & time_value, const float * elem_var_vals, const int & ele_var_index, const int & block_index)
{
	stringstream oss;
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::write_element_variable_to_output()";
	oss << "Writing element variable to output";
	progress_message(&oss,method_name);
#endif
	int error;
	if(!my_output_initialized)
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

template< class field_type >
bool gather_field_data( unsigned expected_num_rel, const field_type & field ,
                        const stk::mesh::Entity & entity ,
                        typename stk::mesh::FieldTraits< field_type >::data_type * dst,
                        stk::mesh::EntityRank entity_rank,
                        const int dim)
{
  typedef typename stk::mesh::FieldTraits< field_type >::data_type T ;

  stk::mesh::PairIterRelation rel = entity.relations( entity_rank );

  bool result = expected_num_rel == (unsigned) rel.size();

  if ( result ) {
    // Iterate over field data for each related entity and copy data
    // into src for one entity at a time
    T * const dst_end = dst + dim * expected_num_rel ;
    for ( ; dst < dst_end ; ++rel , dst += dim ) {
      const T* src = field_data( field , * rel->entity() );
      if (!src) {
        break;
      }
      // FIXME:
      if(dim == 1)
      	stk::Copy<1>( dst , src );
      if(dim == 2)
      	stk::Copy<2>( dst , src );
      if(dim == 3)
      	stk::Copy<3>( dst , src );
    }
    result = dst == dst_end ;
  }
  return result ;
}


bool
Mesh_Manager::verify_coordinates_field(const stk::mesh::STK_Mesh & mesh )
{
	stringstream oss;
#ifdef DEBUG_OUTPUT
	string method_name = "Mesh_Manager::verify_coordinates_field()";
	oss << "Verifying the coordinates in the STK mesh...";
	progress_message(&oss,method_name);
#endif

	bool result = true;

	const stk::mesh::VectorFieldType & coordinates_field = mesh.my_coordinates_field ;
	const stk::mesh::BulkData & bulkData = mesh.my_bulkData ;

	// All element buckets:
	const std::vector<stk::mesh::Bucket*> & elem_buckets =
			bulkData.buckets( stk::mesh::fem::element_rank(mesh.my_fem) );

	// Verify coordinates_field by gathering the nodal coordinates
	// from each element's nodes.
	for ( std::vector<stk::mesh::Bucket*>::const_iterator
			element_bucket_it = elem_buckets.begin();
			element_bucket_it != elem_buckets.end() ; ++element_bucket_it ) {

		const stk::mesh::Bucket& bucket = **element_bucket_it;
		const size_t num_buckets = bucket.size();
		const CellTopologyData * cellTopologyData =
				stk::mesh::fem::get_cell_topology(bucket).getCellTopologyData();

		const int num_nodes = cellTopologyData->node_count;
		const int & dim = mesh.my_spatial_dimension;
		double elem_coord[ num_nodes ][ dim ];

		for( size_t bucket_index = 0; bucket_index < num_buckets; ++bucket_index) {
			const stk::mesh::Entity & elem = bucket[bucket_index] ;

#ifdef DEBUG_OUTPUT
			oss << "Element " << elem.identifier();
			sub_sub_progress_message(&oss);
#endif

			const bool gather_result =
					gather_field_data
					( num_nodes, coordinates_field , elem , & elem_coord[0][0], NODE_RANK, dim );


			if ( gather_result == false ) {
				oss << "verify_coordinates_field() gather was not successful";
				error_message(&oss);
				exit(1);
			}

#ifdef DEBUG_OUTPUT
			for (int node_index=0 ; node_index<num_nodes ; ++node_index )
			{
				log() << "                   node " << node_index + 1 << ": ";
				for (int coord_index=0 ; coord_index<dim ; ++coord_index) {
					log() << "[" << coord_index << "] = " << elem_coord[node_index][coord_index] << " ";
				}
				log() << endl;
			}
#endif
		}
	}
	return result;
}



