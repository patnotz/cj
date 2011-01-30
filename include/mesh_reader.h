/*
 * messages.h
 *
 *  Created on: Jan 30, 2011
 *      Author: dzturne1
 */

#ifndef MESH_READER_H_
#define MESH_READRER_H_

using namespace std;

class Mesh_Reader {
public:
  void read_mesh();
  void set_file_name(char * file_name){my_file_name = file_name;}
private:
  char * my_file_name;
};

#endif /* MESH_READER_H_ */
