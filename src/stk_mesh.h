/*
 * stk_mesh.h
 *
 *  Created on: Feb 6, 2011
 *      Author: dzturne1
 */

#ifndef STK_MESH_H_
#define STK_MESH_H_

#include <vector>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/TopologicalMetaData.hpp>


enum { field_data_chunk_size = 10 };

namespace stk {
namespace mesh {

typedef stk::mesh::Field<double,stk::mesh::Cartesian> VectorFieldType ;
//typedef stk::mesh::Field<double>                      ScalarFieldType ;

class STK_Mesh
{
public:
  ~STK_Mesh();

  STK_Mesh( stk::ParallelMachine comm, const int & spatial_dim);

  const int                      my_spatial_dimension;
  stk::mesh::MetaData            my_metaData;
  stk::mesh::BulkData            my_bulkData;
  stk::mesh::TopologicalMetaData my_topData;
  PartVector my_parts;
  VectorFieldType &              my_coordinates_field;
//  ScalarFieldType &              my_volume_field;
};

} //namespace mesh
} //namespace stk


#endif /* STK_MESH_H_ */
