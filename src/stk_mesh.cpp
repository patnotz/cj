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
{
	// Put the coordinates and temperature field on all nodes
	Part & universal = my_metaData.universal_part();
	put_field(my_coordinates_field,NODE_RANK,universal,my_spatial_dimension);
	// commit will be called by the mesh_manager when it has populated the parts
}

STK_Mesh::~STK_Mesh()
{ }



}
}
