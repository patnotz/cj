#include <stk_mesh.h>
//using stk::mesh::fem::NODE_RANK;

namespace stk {
namespace mesh {

STK_Mesh::STK_Mesh(stk::ParallelMachine comm, const int & spatial_dim) :
    my_spatial_dimension(spatial_dim), my_fem_metaData(my_spatial_dimension), my_bulkData(
      fem::FEMMetaData::get_meta_data(my_fem_metaData), comm), my_elem_rank(
      my_fem_metaData.element_rank()), my_node_rank(
      my_fem_metaData.node_rank()), my_coordinates_field(
      my_fem_metaData.declare_field<VectorFieldType>("coordinates"))
{
  // Put the coordinates and temperature field on all nodes
  Part & universal = my_fem_metaData.universal_part();
  put_field(my_coordinates_field, my_node_rank, universal, my_spatial_dimension);
  // commit will be called by the mesh_manager when it has populated the parts
}

STK_Mesh::~STK_Mesh()
{
}

}
}
