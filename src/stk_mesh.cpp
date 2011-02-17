/*
 * stk_mesh.cpp
 *
 *  Created on: Feb 6, 2011
 *      Author: dzturne1
 */
#include <../include/stk_mesh.h>

namespace stk {
namespace mesh {

STK_Mesh::STK_Mesh( stk::ParallelMachine comm , const int & spatial_dim)
  : my_spatial_dimension(spatial_dim)
  , my_metaData(TopologicalMetaData::entity_rank_names(my_spatial_dimension))
  , my_bulkData(my_metaData,comm,field_data_chunk_size)
  , my_topData(my_metaData,my_spatial_dimension)
  , my_coordinates_field( my_metaData.declare_field< VectorFieldType >("coordinates"))
{
	// Put the coordinates and temperature field on all nodes
	stk::mesh::Part & universal = my_metaData.universal_part();
	stk::mesh::put_field(my_coordinates_field,my_topData.node_rank,universal,my_spatial_dimension);
	// commit will be called by the mesh_manager when it has populated the parts
}

STK_Mesh::~STK_Mesh()
{ }



}
}
