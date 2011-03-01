/*
 * messages.h
 *
 *  Created on: Jan 30, 2011
 *      Author: dzturne1
 */

#ifndef MESSAGES_H_
#define MESSAGES_H_

#include <iostream>
#include <string>

void start_message();
void success_message();
void progress_message(std::stringstream * oss, std::string & name);
void sub_progress_message(std::stringstream * oss);
void sub_sub_progress_message(std::stringstream * oss);
void error_message(std::stringstream * oss);
void printStatus(bool status, std::stringstream * oss);

#endif /* MESSAGES_H_ */
