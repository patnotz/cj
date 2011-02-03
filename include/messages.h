/*
 * messages.h
 *
 *  Created on: Jan 30, 2011
 *      Author: dzturne1
 */

#include <iostream>

#ifndef MESSAGES_H_
#define MESSAGES_H_

using namespace std;

void start_message();
void success_message();
void progress_message(stringstream * oss, string & name);
void sub_progress_message(stringstream * oss);
void sub_sub_progress_message(stringstream * oss);
void error_message(stringstream * oss);

#endif /* MESSAGES_H_ */
