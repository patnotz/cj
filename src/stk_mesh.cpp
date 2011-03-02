/*
 * stk_mesh.cpp
 *
 *  Created on: Feb 6, 2011
 *      Author: dzturne1
 */
#include <stk_mesh.h>
using stk::mesh::fem::NODE_RANK;

namespace stk {
namespace mesh {

STK_Mesh::STK_Mesh( stk::ParallelMachine comm , const int & spatial_dim)
  : my_spatial_dimension(spatial_dim)
  , my_metaData(fem::entity_rank_names(my_spatial_dimension))
  , my_bulkData(my_metaData,comm,field_data_chunk_size)
  , my_fem(my_metaData,my_spatial_dimension)
  , my_elem_rank(fem::element_rank(my_fem))
  , my_coordinates_field( my_metaData.declare_field< VectorFieldType >("coordinates"))
  , my_temperature_field( my_metaData.declare_field< ScalarFieldType >("temperature"))
  , my_volume_field( my_metaData.declare_field< ScalarFieldType >("volume"))
{
	// Put the coordinates and temperature field on all nodes
	stk::mesh::Part & universal = my_metaData.universal_part();
	stk::mesh::put_field(my_coordinates_field,NODE_RANK,universal,my_spatial_dimension);
	stk::mesh::put_field(my_temperature_field,NODE_RANK,universal);
	stk::mesh::put_field(my_volume_field,my_elem_rank,universal);
	// commit will be called by the mesh_manager when it has populated the parts
}

STK_Mesh::~STK_Mesh()
{ }



}
}
