#ifndef STK_MESH_H_
#define STK_MESH_H_

#include <vector>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>

#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

enum
{
  field_data_chunk_size = 10
};

namespace stk {
namespace mesh {

typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorFieldType;
typedef stk::mesh::Field<double> ScalarFieldType;

class STK_Mesh
{
public:
  ~STK_Mesh();

  STK_Mesh(stk::ParallelMachine comm, const int & spatial_dim);

  const int my_spatial_dimension;
  fem::FEMMetaData my_fem_metaData;
  BulkData my_bulkData;
  PartVector my_parts;
  const EntityRank my_elem_rank;
  const EntityRank my_node_rank;

  VectorFieldType & my_coordinates_field;

};

} //namespace mesh
} //namespace stk

#endif /* STK_MESH_H_ */
