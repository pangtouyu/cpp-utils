#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "clock.h"

double Clock::difftime(time_t t1, time_t t2)    // DON'T do this in BC 3.1 
						// Needed only for UNIX boxes
  {
  return t1 - t2;	// NOT t2 - t1
  }

double Clock::elasped_time(void)
  {
  return difftime(time(NULL), start_time);
  }

// given a duration in seconds, prints it as as hh:mm:ss
void Clock::print_time(double time, FILE *fp)
  {
  long hour, min, sec;
  long _sec;

  sec = (long) floor(time + 0.5);
  min = sec/60;
  sec %= 60;
  hour = min/60;
  min %= 60;

  if (fp == NULL)
    printf("%ld:%02ld:%02ld", hour, min, sec);
  else
    fprintf(fp, "%ld:%02ld:%02ld", hour, min, sec);

  return;
  }

void Clock::print_elasped(FILE *fp)
  {
  if (fp == NULL) printf("Elasped time: ");
  else fprintf(fp, "Elasped time: ");
  print_time(elasped_time(), fp);
  if (fp == NULL) printf("\n");
  else fprintf(fp, "\n");
  }

void Clock::print_ETA(float fract_done, FILE *fp)
  {
  double elasped, remaining, tot;

  if (name[0] != 0)
    if (fp == NULL) printf("%s: ", name);
    else fprintf(fp, "%s: ", name);

  elasped = elasped_time();

  if (fract_done > 0)
	  {
	  tot = elasped / fract_done;
	  remaining = tot - elasped;
	  }

  if (fp == NULL) 
    {
    printf("Percent done: %.2f%%    ", 100.0*fract_done);
    printf("Elasped: ");
    print_time(elasped, fp);
    printf("    Remaining: ");
    }
  else
    {
    fprintf(fp, "Percent done: %.2f%%    ", 100.0*fract_done);
    fprintf(fp, "Elasped: ");
    print_time(elasped, fp);
    fprintf(fp, "    Remaining: ");
    }

  if (fract_done > 0)
    print_time(remaining, fp);
  else
    {
    if (fp == NULL) printf("Unknown");
    else fprintf(fp, "Unknown");
    }

  if (fp == NULL) printf("\n");
  else fprintf(fp, "\n");

  return;
  }

void Clock::print_ETA(float amt_done, float tot_work, FILE *fp)
  {
  print_ETA(amt_done/tot_work, fp);
  }

void Clock::reset(void)
  {
  start_time = time(NULL);

  return;
  }

Clock::Clock(const char *in_name)
  {
  name = (char*)malloc(sizeof(char) * (strlen(in_name)+1));
  strcpy(name,in_name);

  start_time = time(NULL);

  return;
  }

Clock::~Clock()
  {
  free(name);

  return;
  }
