#include <iostream>
#include <../include/mesh_reader.h>
#include <../include/messages.h>

//#include "exodusII.h"
//#include "netcdf.h"


using namespace std;

int main()
{

	start_message();

	// This information is hard coded for now, until we have an input reader:
	char* mesh_file_name = "/Users/dzturne1/Documents/dzturne1/Research/cj/problems/unit_2d/mesh.g";

	Mesh_Reader mesh_reader;
	mesh_reader.set_file_name(mesh_file_name);
	mesh_reader.read_mesh();

}


