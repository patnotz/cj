#ifndef ENUMS_H_
#define ENUMS_H_

#include <stdlib.h>
#include <string>
#include <iostream>
//using namespace std;

enum Physics_Type
{
 PERIDYNAMICS = 0,
 FINITE_ELEMENT
};

std::string tostring(const Physics_Type & physics_type);
Physics_Type string_to_physics_type(const std::string & input_string);

#endif /* ENUMS_H_ */
