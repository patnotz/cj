#include <iostream>

#include <drive_simulation.h>
#include <log.h>

int main( int argc, char * argv[] )
{
	Log & log = std::cout;
	const int error_code = drive_simulation(log, argc, argv);
	return error_code;
}
