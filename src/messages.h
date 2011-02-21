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
#include <log.h>

void start_message(Log & log);
void success_message(Log & log);
void progress_message(Log & log, std::stringstream * oss, std::string & name);
void sub_progress_message(Log & log, std::stringstream * oss);
void sub_sub_progress_message(Log & log, std::stringstream * oss);
void error_message(Log & log, std::stringstream * oss);
void printStatus(Log & log, bool status, std::stringstream * oss);

#endif /* MESSAGES_H_ */
