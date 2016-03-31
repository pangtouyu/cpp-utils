#ifndef _READDATA_H
#define _READDATA_H

#include <values.h>

#include "array.h"
#include "matrix.h"

#define INVALID_VALUE (-1)

// Note on variance:
// all variances here are "n" variance, NOT "n-1" variance

class Data
  {
  public:
	int nSamples;
	int nAttributes;		// target is NOT included as an attribute here

	Array<int> attIsContinuous;		// nAttributes long
	int targetIsContinuous;
	Array<int> hasMissingValues;	// nAttributes long

	// the *NValues term will be invalid for continuous attributes

	Array<int> attNValues;		// nAttributes long (note missing values do not 
					//     count as an extra type of value)

	Vector target;			// note target is not allowed to have missing values
	int targetNValues;

	Matrix dataSample;		// nSamples int arrays of length nAttributes each

	Data(int nSamples_, int nAttributes_)
		: nSamples(nSamples_), nAttributes(nAttributes_),
		  attIsContinuous(nAttributes_), hasMissingValues(nAttributes_),
		  attNValues(nAttributes_), target(nSamples_),
		  dataSample(nSamples_, nAttributes_)
		{ }

	// this constructor is meant to be used to creating train/test splits from
	// a Data thing containing all the data. 
	Data(int nSamples_, const Data &allData)
		: nSamples(nSamples_), nAttributes(allData.nAttributes),
		  attIsContinuous(allData.attIsContinuous), 
		  hasMissingValues(allData.hasMissingValues),
		  attNValues(allData.attNValues), 
		  target(nSamples_), dataSample(nSamples_, allData.nAttributes),
		  targetIsContinuous(allData.targetIsContinuous), 
		  targetNValues(allData.targetNValues) 
		{ }

	~Data() { };

	static double missingValueToken(void) { return MAXDOUBLE; }
	void print(FILE *fp, int nSamplesToPrint);

	int numContinuousAtts(void) { int count=0; for (int ctr=0; ctr < nAttributes; ctr++) 
							if (attIsContinuous[ctr]) count++;
				      return count; }
	int inputHasAnyMissingValues(void) { for (int ctr=0; ctr < nAttributes; ctr++)
					    if (hasMissingValues[ctr]) return 1;
					return 0; }

	// calculates some simple statistics on attributes or the target 
	void calculateAttSummary(int att, double *mean, double *variance);
	void calculateTargetSummary(double *mean, double *variance);
	void shuffleSample(void);
  };

#define DATADIR "/opt/data/data/"	// note trailing backslash is required

//-------------------------------------------------------

// discardMissing means any data with 1 or more missing attributes will be
// discarded.  Note that this may make the dataset returned inconsistent
// with what dataSummary() returns.
Data* new_readData(const char dataName[], int discardMissing);

Data* new_recodeMultivariateToBinary(const Data &oldData, int addIntercept);

void dataSummary(const char dataName[], int *nSamples, int *nAttributes);
void new_trainTestSplit(const Data &data,
			Data **trainData_p, Data **testData_p,
			int nTest);
Data* new_addIntercept(const Data &oldData);

#endif		// _READDATA_H

