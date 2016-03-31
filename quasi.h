/* Quasi-Random Number generator, an object oriented
 * implementation in C.
 *
 * modified by Andrew Ng to C++ oop (12/96)
 *
 * returns an n-dimensional vector of values in 0.0..1.0
 * maximum n is currently hardwired to 52
 *
 * See W.H. Press and S.A. Teukolsky, 1989, Quasi- (that is,
 * Sub-) Random Numbers, Computers in Physics V3, No. 6,
 * (Nov/Dec 1989), pp. 76-79
 *
 *
 * rcsid: @(#)quasi.h	1.5 10:15:46 4/18/94   EFC
 *
 */

#ifndef QUASI_RANDOM_H_
#define QUASI_RANDOM_H_ 1.5


struct Quasi
  {
  private:
	int dim;
	unsigned long int index;
	unsigned long int *ix;

  public:
	Quasi(int dimension);
	~Quasi();
	int dimension() {return dim;}

	// get an n-dimensional quasi-random number */
	void QuasiRandomNumber(double x[]);
  };


#endif 		// QUASI_RANDOM_H_ 1.5

