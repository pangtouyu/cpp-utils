#ifndef _CLOCK_H
#define _CLOCK_H

#include <stdlib.h>
#include <time.h>

class Clock
	{
	public:
		double elasped_time(void);
		void print_ETA(float fract_done, FILE *fp = (FILE*)NULL);
		void print_ETA(float amt_done, float tot_work, FILE *fp = (FILE*)NULL);
					// same as print_ETA(amt_done/tot_work);
		void print_elasped(FILE *fp = (FILE*)NULL);
		void reset(void);
		Clock(const char *in_name);
		~Clock(); 

	private:
		char *name;
		time_t start_time;
		void print_time(double time, FILE *fp = (FILE*)NULL);

		double difftime(time_t t1, time_t t2);  
	};

#endif		// _CLOCK_H
