#include <string>
#include <sstream>
#include <time.h>
#include <messages.h>
#include <log.h>

using namespace std;

void
start_message()
{
	log() << endl << endl << "  CODE NAME, please don't call me cj" << endl << endl <<
			"  Authors: PK Notz, DZ Turner" << endl << endl;
	log() << "  Version: 1.0" << endl;
	struct tm *current;
	time_t now;
	time(&now);
	current = localtime(&now);
	int year = current->tm_year + 1900;
	int month = current->tm_mon + 1;
	log() << "  Date:    "
			<< month << "/"
			<<  current->tm_mday << "/"
			<< year << endl;
  	log() << "  Time:    "
  			<< current->tm_hour << ":"
  			<< current->tm_min << ":"
  			<< current->tm_sec << endl;
}

void
success_message()
{
	log() << endl << "  %    Simulation completed successfully    %" << endl;
}


void
progress_message(stringstream * oss, string & name)
{
	string message = oss->str();
	log() << endl << "%  " << name  << endl  << "        " <<  message << endl;
	oss->str(""); oss->clear();
}
void
sub_progress_message(stringstream * oss)
{
	string message = oss->str();
	log() << "             " << message << endl;
	oss->str(""); oss->clear();
}
void
sub_sub_progress_message(stringstream * oss)
{
	string message = oss->str();
	log() << "                   " << message << endl;
	oss->str(""); oss->clear();
}

void
error_message(stringstream * oss)
{
	string message = oss->str();
	error() << "% ERROR % " << message << endl;
	oss->str(""); oss->clear();
}

void
printStatus(bool status, stringstream * oss)
{
	string message = oss->str();
  if (status) {
    log() << "  " << message << " [ PASSED ]" << endl;
  }
  else {
    log() << "  " << message << "[ FAILED ]" << endl;
  }
	oss->str(""); oss->clear();
}

