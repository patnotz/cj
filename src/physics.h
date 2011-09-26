#ifndef PHYSICS_H_
#define PHYSICS_H_

#include <Teuchos_RCP.hpp>
#include <enums.h>
#include <stk_mesh.h>
#include <analysis_model.h>

class Physics
{
public:
  Physics(const std::string & analysis_model_str, const Physics_Type physics_type, const std::string & name);
  virtual ~Physics(){};

  const Physics_Type physics_type(){return my_physics_type;}
  const std::string analysis_model_str(){return my_analysis_model_str;}
  const std::string name(){return my_name;}
  Teuchos::RCP<stk::mesh::STK_Mesh> & stk_mesh_ptr(){return my_stk_mesh;}

  virtual void initialize_fields(Analysis_Model & analysis_model)=0;
  void initialize_mesh_parts_and_commit(const Analysis_Model & analysis_model);
  stk::mesh::Part * const
  part_pointer(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh, const std::string & elem_type,
    const std::string & name);
  void print_field_info(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh);
  void populate_mesh_elements(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh, const Analysis_Model & analysis_model);
  void populate_mesh_coordinates(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh,const Analysis_Model & analysis_model);
  bool verify_coordinates_field(Teuchos::RCP<stk::mesh::STK_Mesh> mesh, const Analysis_Model & analysis_model);

  void populate_bogus_scalar_field(Teuchos::RCP<stk::mesh::STK_Mesh>  const mesh,
    stk::mesh::ScalarFieldType & field, const stk::mesh::EntityRank & rank);
  void populate_bogus_vector_field(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh,
    stk::mesh::VectorFieldType & field, const stk::mesh::EntityRank & rank);


protected:
  Physics(const Physics&);
  Physics& operator=(const Physics&);

  std::string my_analysis_model_str;
  Physics_Type my_physics_type;
  std::string my_name;
  Teuchos::RCP<stk::mesh::STK_Mesh> my_stk_mesh;

};

std::vector<std::string> get_field_names(Teuchos::RCP<stk::mesh::STK_Mesh> const mesh,
  stk::mesh::EntityRank entity_rank, unsigned field_rank);


class Peridynamics_Physics : public Physics
{
public:
  Peridynamics_Physics(const std::string & analysis_model_str, const Physics_Type physics_type, const std::string & name);
  virtual ~Peridynamics_Physics(){};

  virtual void initialize_fields(Analysis_Model & analysis_model);

protected:
  Peridynamics_Physics(const Peridynamics_Physics&);
  Peridynamics_Physics& operator=(const Peridynamics_Physics&);

};


class Physics_Factory
{
public:
  Physics_Factory();
  virtual ~Physics_Factory(){}
  virtual Teuchos::RCP<Physics> create(const std::string & analysis_model_str,const Physics_Type physics_type,const std::string & name);

private:
  Physics_Factory(const Physics_Factory&);
  Physics_Factory& operator=(const Physics_Factory&);

protected:

};

#endif /* PHYSICS_H_ */
