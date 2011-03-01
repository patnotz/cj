#include <log.h>
#include <iostream>

std::ostream * Log::p_log = & std::cout;
std::ostream * Log::p_error = & std::cerr;
std::ostream * Log::p_warning = & std::cout;

void Log::setLog(std::ostream & log) {
	p_log = &log;
}

void Log::setError(std::ostream & error) {
	p_error = &error;
}

void Log::setWarning(std::ostream & warning) {
	p_warning = &warning;
}

std::ostream & Log::log() {
	return *p_log;
}

std::ostream & Log::error() {
	return *p_error;
}

std::ostream & Log::warning() {
	return *p_warning;
}

std::ostream & log() {
	return Log::log();
}

std::ostream & error() {
	return Log::error();
}

std::ostream & warning() {
	return Log::warning();
}
