/*
 * messages.cpp
 *
 *  Created on: Jan 30, 2011
 *      Author: dzturne1
 */
#include <string>
#include <sstream>
#include <time.h>
#include <../include/messages.h>
//#include <../include/trace.h>

using namespace std;

void
start_message()
{
	cout << endl << endl << "  CODE NAME, please don't call me cj" << endl << endl <<
			"  Authors: PK Notz, DZ Turner" << endl << endl;
	cout << "  Version: 1.0" << endl;
	struct tm *current;
	time_t now;
	time(&now);
	current = localtime(&now);
	int year = current->tm_year + 1900;
	int month = current->tm_mon + 1;
	cout << "  Date:    " << month << "/" <<  current->tm_mday << "/" << year << endl;
	printf("  Time:    %i:%i:%i\n", current->tm_hour, current->tm_min, current->tm_sec);
}

void
success_message()
{
	cout << endl << "  %    Simulation completed successfully    %" << endl;
}


void
progress_message(stringstream * oss, string & name)
{
	string message = oss->str();
	cout << endl << "%  " << name  << endl  << "        " <<  message << endl;
	oss->str(""); oss->clear();
}
void
sub_progress_message(stringstream * oss)
{
	string message = oss->str();
	cout << "             " << message << endl;
	oss->str(""); oss->clear();
}
void
sub_sub_progress_message(stringstream * oss)
{
	string message = oss->str();
	cout << "                   " << message << endl;
	oss->str(""); oss->clear();
}

void
error_message(stringstream * oss)
{
	string message = oss->str();
	cout << "% ERROR % " << message << endl;
	oss->str(""); oss->clear();
}

