#ifndef _NORMAL_H
#define _NORMAL_H

double normal_dev(double mean, double variance);
double normal_dev_stddev(double mean, double stddev);
double gasdev(void);
inline double standardNormal_dev(void) {return gasdev();}

#endif		/* _NORMAL_H */

