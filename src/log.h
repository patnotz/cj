#ifndef LOG_H_
#define LOG_H_

#include <ostream>

class Log {
public:
	static void setLog(std::ostream & log);
	static void setError(std::ostream & error);
	static void setWarning(std::ostream & warning);
	static std::ostream & log();
	static std::ostream & error();
	static std::ostream & warning();
private:
	static std::ostream * p_log;
	static std::ostream * p_error;
	static std::ostream * p_warning;
};

std::ostream & log();
std::ostream & error();
std::ostream & warning();

//#define DEBUG_OUTPUT

#endif /* LOG_H_ */
