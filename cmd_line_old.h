#ifndef _CMD_LINE_H
#define _CMD_LINE_H
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#define MAX_FLAG_SIZE 64

class cmd_line
  {
  private:
        int argc;
        char **argv;
	int *arg_used;
	char *cmd_line_string;

  public:
		int flag_position(const char* flag) const;
		int get_flag_pos(const char* flag) const {return flag_position(flag);}
		int exist_flag(const char *in_flag) const;
		float get_flag_arg(const char *in_flag, float dfault, int offset=0);
		double get_flag_arg(const char *in_flag, double dfault, int offset=0);
		long get_flag_arg(const char *in_flag, long dfault, int offset=0);
		const char *get_flag_arg(const char *in_flag, const char *dfault, int offset=0);
		const char *cmd_line_str(void) {return cmd_line_string;}
		void assert_all_parms_used(void) const;
		int get_flag_parm(const char *in_flag, int nterms, ...) const;

		cmd_line(int argc, char *argv[]);
		cmd_line::~cmd_line(void);
  };

void start_random_number_from_cline(cmd_line &c_line,
                                int *seed1p=(int*)NULL, int *seed2p=(int*)NULL);

/*
BUGS: 

1. get_flag_parm() doesn't work

*/

#endif		// _CMD_LINE_H
