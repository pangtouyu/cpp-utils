/**************************************

  cmd_line.C - does stuff with
    the command line for you.

  Copyright(c) 1994, Andrew Y. Ng
  All Rights reserved.

**************************************/

#include "cmd_line.h"
#include "Random_Number.h"

cmd_line::cmd_line(int in_argc, char *in_argv[])
  {
  int tot_len;
  int ctr, dex;

  argv=(char**)malloc(sizeof(char*)*(in_argc+1));
  arg_used = (int*)calloc(sizeof(int), in_argc);
  argc=in_argc;
 
  tot_len=0;
  for (ctr=0; ctr < in_argc; ctr++)
	{
	argv[ctr] = (char*)malloc(sizeof(char) * 
			(strlen(in_argv[ctr])+1));
	strcpy(argv[ctr], in_argv[ctr]);

	tot_len += strlen(in_argv[ctr])+1;
	}

  dex=0;
  cmd_line_string = (char*)malloc(sizeof(char) * tot_len);
  for (ctr=0; ctr < in_argc; ctr++)
	{
	strcpy(cmd_line_string+dex, in_argv[ctr]);
	dex += strlen(in_argv[ctr])+1;
	if (ctr != in_argc-1)
		cmd_line_string[dex-1] = ' ';
	}

  return;
  }

cmd_line::~cmd_line(void)
  {
  int ctr;

  for (ctr=0; ctr < argc; ctr++)
	free(argv[ctr]);
  free(arg_used);
  free(argv);
  
  free(cmd_line_string);

  return;
  }

int cmd_line::flag_position(const char *in_flag) const
  {
  char flag[MAX_FLAG_SIZE] = {'-', 0};
  int ctr;

  strcat(flag, in_flag);

  for (ctr=1; ctr < argc; ctr++)
	if (!strcmp(argv[ctr], flag))
		{
		arg_used[ctr] = 1;
		return ctr;
		}

  return 0;
  }

int cmd_line::exist_flag(const char *in_flag) const
  {
  return (flag_position(in_flag) != 0);
  }

float cmd_line::get_flag_arg(const char *in_flag, float dfault, int offset)
// dfault is "default"
  {
  int ctr;
  int flag_pos;
  float temp;

  flag_pos = flag_position(in_flag);
  if (flag_pos == 0)
	return dfault;
  for (ctr=flag_pos+1; ctr <= flag_pos+1+offset; ctr++)
	if (argv[ctr] == NULL)
		return dfault;

  if (sscanf(argv[flag_pos+1+offset], "%f", &temp) == 0)
	return dfault;

  arg_used[flag_pos+1+offset] = 1;
  return temp;
  }

double cmd_line::get_flag_arg(const char *in_flag, double dfault, int offset)
// dfault is "default"
  {
  int ctr;
  int flag_pos;
  double temp;

  flag_pos = flag_position(in_flag);
  if (flag_pos == 0)
	return dfault;
  for (ctr=flag_pos+1; ctr <= flag_pos+1+offset; ctr++)
	if (argv[ctr] == NULL)
		return dfault;

  if (sscanf(argv[flag_pos+1+offset], "%lf", &temp) == 0)
	return dfault;

  arg_used[flag_pos+1+offset] = 1;
  return temp;
  }

long cmd_line::get_flag_arg(const char *in_flag, long dfault, int offset)
// dfault is "default"
  {
  int ctr;
  int flag_pos;
  long temp;

  flag_pos = flag_position(in_flag);
  if (flag_pos == 0)
	return dfault;
  for (ctr=flag_pos+1; ctr <= flag_pos+1+offset; ctr++)
	if (argv[ctr] == NULL)
		return dfault;

  if (sscanf(argv[flag_pos+1+offset], "%ld", &temp) == 0)
		return dfault;

  arg_used[flag_pos+1+offset] = 1;
  return temp;
  }

const char* cmd_line::get_flag_arg(const char *in_flag, const char *dfault, int offset)
// dfault is "default"
  {
  int ctr;
  int flag_pos;
  const char* temp;

  flag_pos = flag_position(in_flag);
  if (flag_pos == 0)
	return dfault;
  for (ctr=flag_pos+1; ctr <= flag_pos+1+offset; ctr++)
	if (argv[ctr] == NULL)
		return dfault;

  temp = argv[flag_pos+1+offset];

  arg_used[flag_pos+1+offset] = 1;
  return temp;
  }

void cmd_line::assert_all_parms_used(void) const
  {
  int ctr;

  // ctr starts at 1 to ignore the command/program name
  for (ctr=1; ctr < argc; ctr++)
	if (!arg_used[ctr])
		{
		fprintf(stderr,"ERROR (cmd_line): Extra argument: %s\n", 
					argv[ctr]);
		exit(-1);
		}

  return;
  }

// called with something like get_flag_parm("time", 2, "%d", &tw, "%d", &tb);
// this function doesn't work. (yet!)
int cmd_line::get_flag_parm(const char *in_flag, int nterms, ...) const
  {
  char flag[MAX_FLAG_SIZE] = {'-', 0};
  int ctr, flag_pos, n_terms_read;
//  int temp_nterms;
  const char* format;
  void *temp_addr;
  va_list ap;

//******************************************************** 
  fprintf(stderr,"sorry, get_flag_parm() doesn't work.\n");
  exit(-1);

  strcat(flag, in_flag);
  if ((flag_pos=flag_position(in_flag))== 0)
	return 0;		// no such flag

//  temp_nterms = 2*nterms+2;	// va_start is a macro, and
				// will complain if I don't do
				// this

// The next line is supposed to be there, but radish's compiler
// doesn't like it somehow, so I don't use it.
//  va_start(ap, (2*nterms+2));

  n_terms_read=0;

  for (ctr=1; ctr <= nterms; ctr++)
	{
	format=va_arg(ap, char*);
	printf("Format %d: %s\n", ctr, format);
	temp_addr=va_arg(ap, void*);
	n_terms_read += sscanf(argv[flag_pos+ctr], format, temp_addr);
	}

  va_end(ap);

  return n_terms_read;
  }

//-------------------------------------------------------

void start_random_number_from_cline(cmd_line &c_line,
                                int *seed1p=(int*)NULL, int *seed2p=(int*)NULL)
  {
  const long DEFAULT_SEED1 = 98239,
             DEFAULT_SEED2 = 19281;

  int seed1, seed2;

  seed1 = c_line.get_flag_arg("seeds", DEFAULT_SEED1, 0);
  seed2 = c_line.get_flag_arg("seeds", DEFAULT_SEED2, 1);

  start_random_number(seed1, seed2);

  if (seed1p != NULL)
        *seed1p = seed1;
  if (seed2p != NULL)
        *seed2p = seed2;

  return;
  }


