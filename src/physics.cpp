#include <physics.h>
#include <enums.h>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>
#include <stk_mesh/fem/Stencils.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <analysis_model.h>
#include <map>
#include <Teuchos_RCP.hpp>
#include <log.h>
#include <analysis_model.h>
#include <messages.h>

using namespace std;

Physics_Factory::Physics_Factory(){}

Teuchos::RCP<Physics> Physics_Factory::create(const std::string & analysis_model_str,const Physics_Type physics_type,const std::string & name)
{
  using Teuchos::rcp;

  if(physics_type==PERIDYNAMICS)
    return rcp(new Peridynamics_Physics(analysis_model_str, physics_type,name));
  else if (physics_type==FINITE_ELEMENT)
    return rcp(new Peridynamics_Physics(analysis_model_str, physics_type,name));
  else
  {
    cerr << "Error: unknown physics type cannot be created: " << tostring(physics_type) << std::endl;
    exit(1);
  }
}

Physics::Physics(const std::string & analysis_model_str,const Physics_Type physics_type,const std::string & name) :
    my_analysis_model_str(analysis_model_str),
    my_physics_type(physics_type),
    my_name(name)
{
}

void Physics::initialize_mesh_parts_and_commit(const Analysis_Model & analysis_model)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Physics::initialize_mesh_parts()";
  oss << "Initializing mesh parts";
  progress_message(&oss, method_name);
#endif
  if (!analysis_model.input_initialized())
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
  for (int i = 0; i < analysis_model.num_blocks(); ++i)
  {
    const string elem_type = analysis_model.elem_types().find(analysis_model.block_ids()[i])->second;
    oss << "block_" << analysis_model.block_ids()[i] << "_" << elem_type;
    const string part_name = oss.str();
    oss.str("");
    oss.clear();

    stk::mesh::Part * part_ptr = part_pointer(my_stk_mesh, elem_type, part_name);
    my_stk_mesh->my_parts.push_back(part_ptr);
  }
#ifdef DEBUG_OUTPUT
  for (int i = 0; i < analysis_model.num_blocks(); ++i)
  {
    oss << my_stk_mesh->my_parts[i]->name();
    sub_sub_progress_message(&oss);
  }
#endif
  my_stk_mesh->my_fem_metaData.commit(); // this has to be called before we can modify the mesh
#ifdef DEBUG_OUTPUT
  oss << "Mesh committed, adding the elements";
  sub_progress_message(&oss);
#endif
}

std::vector<std::string> get_field_names(
  Teuchos::RCP<stk::mesh::STK_Mesh> const mesh, stk::mesh::EntityRank entity_rank,
  unsigned field_rank)
{
  std::vector<std::string> nodal_vars_vec;
  const stk::mesh::FieldVector & fields = mesh->my_fem_metaData.get_fields();
  unsigned nfields = fields.size();
  //std::cout << "Number of fields = " << fields.size() << std::endl;
  for (unsigned ifld = 0; ifld < nfields; ifld++)
  {
    stk::mesh::FieldBase *field = fields[ifld];
    //std::cout << "Field[" << ifld << "]= " << field->name() << " rank= " << field->rank() << std::endl;
    //std::cout << "info>    " << *field << std::endl;
    if (field->rank() == field_rank) // scalar or vector
    {
      unsigned nfr = field->restrictions().size();
      //std::cout << "number of field restrictions= " << nfr << std::endl;
      for (unsigned ifr = 0; ifr < nfr; ifr++)
      {
        const stk::mesh::FieldRestriction & fr = field->restrictions()[ifr];
        //std::cout << fr.key.rank();
        if (fr.entity_rank() == entity_rank) // nodal, element, etc.
        {
          string field_name = field->name();
          //std::cout << "this is the field name: " << field_name << endl;
          nodal_vars_vec.push_back(field_name);
          //std::cout << "Field[" << ifld << "]= " << field->name() << std::endl;
        }
      }
    }
  }
  return nodal_vars_vec;
}


void Physics::print_field_info(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh)
{
  log() << endl
        << "  --------------------------------------------------------------------------------------"  << endl
        << "    Physics: " << my_name << endl
        << "    Physics type: " << tostring(my_physics_type) << endl
        << "    Analysis model: " << my_analysis_model_str << endl;

  log() << "  =============== REGISTERED NODAL SCALAR FIELDS ==============="
      << endl;
  std::vector<std::string> nodal_scalar_fields = get_field_names(mesh,
    mesh->my_node_rank, 0);
  for (int i = 0; i < nodal_scalar_fields.size(); ++i)
  {
    log() << "  " << nodal_scalar_fields[i] << endl;
  }
  log() << "  =============== REGISTERED NODAL VECTOR FIELDS ==============="
      << endl;
  std::vector<std::string> nodal_vector_fields = get_field_names(mesh,
    mesh->my_node_rank, 1);
  for (int i = 0; i < nodal_vector_fields.size(); ++i)
  {
    log() << "  " << nodal_vector_fields[i] + "_x" << endl;
    log() << "  " << nodal_vector_fields[i] + "_y" << endl;
    if (mesh->my_spatial_dimension > 2) log() << "  "
        << nodal_vector_fields[i] + "_z" << endl;
  }
  log() << "  ============== REGISTERED ELEMENT SCALAR FIELDS =============="
      << endl;
  std::vector<std::string> element_scalar_fields = get_field_names(mesh,
    mesh->my_elem_rank, 0);
  for (int i = 0; i < element_scalar_fields.size(); ++i)
  {
    log() << "  " << element_scalar_fields[i] << endl;
  }
  log() << "  ============== REGISTERED ELEMENT VECTOR FIELDS =============="
      << endl;
  std::vector<std::string> element_vector_fields = get_field_names(mesh,
    mesh->my_elem_rank, 1);
  for (int i = 0; i < element_vector_fields.size(); ++i)
  {
    log() << "  " << element_vector_fields[i] << endl;
  }
  log() << "  =============================================================="
      << endl;

  log() <<  "  --------------------------------------------------------------------------------------" << endl;

}

void Physics::populate_mesh_elements(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh, const Analysis_Model & analysis_model)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Physics::populate_mesh_elements()";
  oss << "Populating mesh elements";
  progress_message(&oss, method_name);
#endif
  mesh->my_bulkData.modification_begin(); // Begin modifying the mesh
  int ele_map_index = 0;
  for (int block = 0; block < analysis_model.num_blocks(); ++block)
  {
    int * connectivity = analysis_model.connectivities().find(analysis_model.block_ids()[block])->second;
    const string elem_type = analysis_model.elem_types().find(analysis_model.block_ids()[block])->second;
    for (int ele = 0; ele < analysis_model.num_elem_in_block(block); ++ele)
    {
      stk::mesh::EntityId node_ids[analysis_model.num_nodes_per_elem()[block]];
      stk::mesh::EntityId elem_id = analysis_model.elem_map()[ele_map_index];

      // Note declare_element expects a cell topology
      // to have been attached to the part
      analysis_model.map_node_ids(block, ele, node_ids, elem_type, connectivity);
#ifdef DEBUG_OUTPUT
      oss << "ExodusII: ";
      int base = 0;
      if (elem_type == "QUAD4")
      {
        base = (ele) * 4;
      }
      else if (elem_type == "HEX8")
      {
        base = (ele) * 8;
      }
      else if (elem_type == "TETRA4" || elem_type == "TETRA")
      {
        base = (ele) * 4;
      }
      else if (elem_type == "TRI3")
      {
        base = (ele) * 3;
      }
      for (int node = 0; node < analysis_model.num_nodes_per_elem()[block]; ++node)
      {
        oss << connectivity[base + node] << ",";
      }
      oss << " -> STK mesh: ";
      for (int node = 0; node < analysis_model.num_nodes_per_elem()[block]; ++node)
      {
        oss << node_ids[node] << ",";
      }
      sub_sub_progress_message(&oss);
#endif
#ifdef DEBUG_OUTPUT
      oss << "Adding element id: " << elem_id << " in block "
          << mesh->my_parts[block]->name();
      sub_sub_progress_message(&oss);
#endif
      stk::mesh::fem::declare_element(mesh->my_bulkData, *mesh->my_parts[block],
        elem_id, node_ids);
      ele_map_index++;
    }
  }
  // Done modifying the mesh.
  // Modifications on the local parallel process are communicated
  // among processes, verified for consistency, and changes to
  // parallel shared/ghosted mesh entities are synchronized.
  mesh->my_bulkData.modification_end();
}

void Physics::populate_mesh_coordinates(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh, const Analysis_Model & analysis_model)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Physics::populate_mesh_coorindates()";
  oss << "Populating STK mesh";
  progress_message(&oss, method_name);
#endif
  const std::vector<stk::mesh::Bucket*> & node_buckets =
      mesh->my_bulkData.buckets(mesh->my_node_rank);

  for (std::vector<stk::mesh::Bucket*>::const_iterator node_bucket_it =
      node_buckets.begin(); node_bucket_it != node_buckets.end();
      ++node_bucket_it)
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
    stk::mesh::BucketArray<stk::mesh::VectorFieldType> coordinates_array(
      mesh->my_coordinates_field, bucket);
    const int num_sets_of_coords = coordinates_array.dimension(1); //this is linked to the bucket size
    for (int i = 0; i < num_sets_of_coords; ++i)
    {
      const unsigned node_id = bucket[i].identifier();
      analysis_model.map_node_coordinates(node_id, &coordinates_array(0, i));
    }
  }
  //Now that x, y, and z are in the stk mesh field we can delete the temp data storage
  delete analysis_model.x_coord();
  delete analysis_model.y_coord();
  if (analysis_model.num_dim() > 2) delete analysis_model.z_coord();
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


bool Physics::verify_coordinates_field(Teuchos::RCP<stk::mesh::STK_Mesh> mesh, const Analysis_Model & analysis_model)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Physics::verify_coordinates_field()";
  oss << "Verifying the coordinates in the STK mesh...";
  progress_message(&oss, method_name);
#endif

  bool result = true;

  const stk::mesh::VectorFieldType & coordinates_field =
      mesh->my_coordinates_field;
  const stk::mesh::BulkData & bulkData = mesh->my_bulkData;

  // All element buckets:
  const std::vector<stk::mesh::Bucket*> & elem_buckets = bulkData.buckets(
    mesh->my_elem_rank);

  // Verify coordinates_field by gathering the nodal coordinates
  // from each element's nodes.
  for (std::vector<stk::mesh::Bucket*>::const_iterator element_bucket_it =
      elem_buckets.begin(); element_bucket_it != elem_buckets.end();
      ++element_bucket_it)
      {

    const stk::mesh::Bucket& bucket = **element_bucket_it;
    const size_t num_buckets = bucket.size();
    const CellTopologyData * cellTopologyData =
        stk::mesh::fem::get_cell_topology(bucket).getCellTopologyData();

    const int num_nodes = cellTopologyData->node_count;
    const int & dim = mesh->my_spatial_dimension;
    double elem_coord[num_nodes][dim];

    for (size_t bucket_index = 0; bucket_index < num_buckets; ++bucket_index)
    {
      const stk::mesh::Entity & elem = bucket[bucket_index];

#ifdef DEBUG_OUTPUT
      oss << "Element " << elem.identifier();
      sub_sub_progress_message(&oss);
#endif

      const bool gather_result = gather_field_data(num_nodes, coordinates_field,
        elem, &elem_coord[0][0], mesh->my_node_rank, dim);

      if (gather_result == false)
      {
        oss << "verify_coordinates_field() gather was not successful";
        error_message(&oss);
        exit(1);
      }

#ifdef DEBUG_OUTPUT
      for (int node_index = 0; node_index < num_nodes; ++node_index)
      {
        log() << "                   node " << node_index + 1 << ": ";
        for (int coord_index = 0; coord_index < dim; ++coord_index)
        {
          log() << "[" << coord_index << "] = "
              << elem_coord[node_index][coord_index] << " ";
        }
        log() << endl;
      }
#endif
    }
  }
  return result;
}

void Physics::populate_bogus_scalar_field(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh,
  stk::mesh::ScalarFieldType & field, const stk::mesh::EntityRank & rank)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Physics::populate_bogus_scalar_field()";
  oss << "Populating bogus scalar field " << field.name();
  progress_message(&oss, method_name);
#endif
  const std::vector<stk::mesh::Bucket*> & buckets = mesh->my_bulkData.buckets(
    rank);

  for (std::vector<stk::mesh::Bucket*>::const_iterator bucket_it =
      buckets.begin(); bucket_it != buckets.end(); ++bucket_it)
  {
    const stk::mesh::Bucket & bucket = **bucket_it;
    stk::mesh::BucketArray<stk::mesh::ScalarFieldType> field_array(field,
      bucket);
    const int num_items_in_bucket = field_array.dimension(0); //this is linked to the bucket size
    for (int i = 0; i < num_items_in_bucket; ++i)
    {
      const unsigned item_id = bucket[i].identifier();
      if (rank == mesh->my_elem_rank) field_array(i) = item_id * 1150.6;
      else
        field_array(i) = item_id * 4.6890;
    }
  }
}

void Physics::populate_bogus_vector_field(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh,
  stk::mesh::VectorFieldType & field, const stk::mesh::EntityRank & rank)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Physics::populate_bogus_vector_field()";
  oss << "Populating bogus vector field " << field.name();
  progress_message(&oss, method_name);
#endif
  const std::vector<stk::mesh::Bucket*> & buckets = mesh->my_bulkData.buckets(
    rank);

  for (std::vector<stk::mesh::Bucket*>::const_iterator bucket_it =
      buckets.begin(); bucket_it != buckets.end(); ++bucket_it)
  {
    const stk::mesh::Bucket & bucket = **bucket_it;
    stk::mesh::BucketArray<stk::mesh::VectorFieldType> field_array(field,
      bucket);
    const int num_components = field_array.dimension(0); // this should be linked to the spatial dim
    const int num_items_in_bucket = field_array.dimension(1); //this is linked to the bucket size
    for (int i = 0; i < num_items_in_bucket; ++i)
    {
      const unsigned id = bucket[i].identifier();
      double * values = &field_array(0, i);
      values[0] = id * 100000.0;
      values[1] = id * 100100.0;
      if (num_components > 2) values[2] = id * 100200.0;
    }
  }
}



Peridynamics_Physics::Peridynamics_Physics(const std::string & analysis_model_str,const Physics_Type physics_type,const std::string & name) :
    Physics(analysis_model_str,physics_type,name)
{
}

void
Peridynamics_Physics::initialize_fields(Analysis_Model & analysis_model)
{
  stringstream oss;
#ifdef DEBUG_OUTPUT
  string method_name = "Peridynamics_Physics::initialize_fields()";
  oss << "Initializing fields for physics: " << my_name;
  progress_message(&oss, method_name);
#endif

  stk::mesh::Part & universal = my_stk_mesh->my_fem_metaData.universal_part();

#ifdef DEBUG_OUTPUT
  oss << "Got a universal part from the stk mesh";
  sub_progress_message(&oss);
#endif


  //    stk::mesh::ScalarFieldType & temperature_field =
  //        stk_mesh.my_fem_metaData.declare_field<stk::mesh::ScalarFieldType>(
  //        "temperature");
  //    stk::mesh::put_field(temperature_field, stk_mesh.my_node_rank, universal);
  //
  //    stk::mesh::VectorFieldType & velocity_field =
  //        stk_mesh.my_fem_metaData.declare_field<stk::mesh::VectorFieldType>(
  //          "velocity");
  //    stk::mesh::put_field(velocity_field, stk_mesh.my_node_rank, universal,
  //      stk_mesh.my_spatial_dimension);
  //
  //    stk::mesh::ScalarFieldType & volume_field =
  //        stk_mesh.my_fem_metaData.declare_field<stk::mesh::ScalarFieldType>(
  //          "volume"); //
  //    stk::mesh::put_field(volume_field, stk_mesh.my_elem_rank, universal);
  //
  //    stk::mesh::VectorFieldType & acceleration_field =
  //        stk_mesh.my_fem_metaData.declare_field<stk::mesh::VectorFieldType>(
  //          "acceleration");
  //    stk::mesh::put_field(acceleration_field, stk_mesh.my_elem_rank, universal);


  stk::mesh::ScalarFieldType & volume_field =
      my_stk_mesh->my_fem_metaData.declare_field<stk::mesh::ScalarFieldType>(
        "volume"); //
  stk::mesh::put_field(volume_field, my_stk_mesh->my_elem_rank, universal);

#ifdef DEBUG_OUTPUT
  oss << "Put field onto the stk mesh successful: volume";
  sub_progress_message(&oss);
#endif


  initialize_mesh_parts_and_commit(analysis_model);
  print_field_info(my_stk_mesh);

  // put the elements and connectivity in the stk_mesh
  populate_mesh_elements(my_stk_mesh,analysis_model);
  // put the coordinates as a field on the stk_mesh
  populate_mesh_coordinates(my_stk_mesh,analysis_model);
  bool local_status = true;
  local_status = verify_coordinates_field(my_stk_mesh,analysis_model);

  oss << "Verifying the STK mesh coordinates field: ";
  printStatus(local_status, &oss);

  //    pos->second->populate_bogus_scalar_field(&stk_mesh, temperature_field,
  //      stk_mesh.my_node_rank);
  //    // volume
  //    pos->second->populate_bogus_scalar_field(&stk_mesh, volume_field,
  //      stk_mesh.my_elem_rank);
  //    // velocity
  //    pos->second->populate_bogus_vector_field(&stk_mesh, velocity_field,
  //      stk_mesh.my_node_rank);
  //    // acceleration
  //    pos->second->populate_bogus_vector_field(&stk_mesh, acceleration_field,
  //      stk_mesh.my_elem_rank);


  // POPULATE BOGUS FIELDS FOR NOW
  //  volume
  populate_bogus_scalar_field(my_stk_mesh, volume_field,
    my_stk_mesh->my_elem_rank);

  const char * title = my_name.c_str();
  analysis_model.initialize_output(title, my_stk_mesh);

  int time_step = 1;
  int num_time_steps = 1;
  for (int i = 0; i < num_time_steps; i++)
  {
    float time_value = (float) (i + 1) / 100.0;
    analysis_model.write_time_step_info(time_step, time_value);
    analysis_model.write_element_scalar(time_step, time_value, my_stk_mesh,
      volume_field);
    analysis_model.update_output();
    time_step++;
  }
}

stk::mesh::Part * const
Physics::part_pointer(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh,
  const string & elem_type, const string & name)
{
  stringstream oss;
  stk::mesh::Part * part_ptr;

  if (elem_type == "QUAD4")
  {
    stk::mesh::Part & part = stk::mesh::fem::declare_part<
        shards::Quadrilateral<4> >(mesh->my_fem_metaData, name);
    part_ptr = &part;
  }
  else if (elem_type == "HEX8")
  {
    stk::mesh::Part & part =
        stk::mesh::fem::declare_part<shards::Hexahedron<8> >(
          mesh->my_fem_metaData, name);
    part_ptr = &part;
  }
  else if (elem_type == "TETRA4" || elem_type == "TETRA")
  {
    stk::mesh::Part & part =
        stk::mesh::fem::declare_part<shards::Tetrahedron<4> >(
          mesh->my_fem_metaData, name);
    part_ptr = &part;
  }
  else if (elem_type == "TRI3")
  {
    stk::mesh::Part & part = stk::mesh::fem::declare_part<shards::Triangle<3> >(
      mesh->my_fem_metaData, name);
    part_ptr = &part;
  }
  else if (elem_type == "SPHERE")
  {
    stk::mesh::Part & part = stk::mesh::fem::declare_part<shards::Particle>(
      mesh->my_fem_metaData, name);
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

