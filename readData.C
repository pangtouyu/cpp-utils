#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <values.h>

#include "array.h"
#include "readData.h"
#include "matrix.h"
#include "misc.h"
#include "normal.h"

//-------------------------------------------------------
// class Data

void Data::calculateAttSummary(int att, double *mean, double *variance)
  {
  double sum = 0, exSqr=0;
  int denom = 0;

  for (int s=0; s < nSamples; s++)
	{
	double val = dataSample[s][att];
	if (val != missingValueToken())
	    {
	    sum += val; 
	    exSqr += sqr(val);
	    denom++;
	    }
	}

  if (denom >= 2)
	{
	*mean = sum/denom;
	exSqr /= denom;
//	double temp = (exSqr - sqr(*mean)) * (double)denom / (denom-1);
	double temp = (exSqr - sqr(*mean));
	*variance = max(temp, 0.0);	// in case roundoff makes "Variance < 0"
	}
  else if (denom == 1)
	{ *mean = sum/denom; *variance = 0.0; }
  else 
	{ *mean = *variance = 0.0; }

  return;
  }

void Data::calculateTargetSummary(double *mean, double *variance)
  {
  double sum = 0, exSqr=0;
  int denom = 0;

  for (int s=0; s < nSamples; s++)
	{
	double val = target[s];
	assert(val != missingValueToken());
	sum += val; 
	exSqr += sqr(val);
	denom++;
	}

  if (denom >= 2)
	{
	*mean = sum/denom;
	exSqr /= denom;
//	double temp = (exSqr - sqr(*mean)) * (double)denom / (denom-1);
	double temp = (exSqr - sqr(*mean));
	*variance = max(temp, 0.0);	// in case roundoff makes "Variance < 0"
	}
  else if (denom == 1)
	{ *mean = sum/denom; *variance = 0.0; }
  else 
	{ *mean = *variance = 0.0; }

  return;
  }

void Data::shuffleSample(void)
  {
  for (int s1=0; s1 < nSamples; s1++)
    {
    int s2 = rand_int(0, nSamples-1);

    // will swap samples s1 and s1;
    for (int att=0; att < nAttributes; att++)
	swap(dataSample[s1][att], dataSample[s2][att]);
    swap(target[s1], target[s2]);
    }

  return;
  }

void Data::print(FILE *fp, int nSamplesToPrint)
  {
  printf("%d samples. %d attributes.\n", nSamples, nAttributes);
  printf("continuous/discrete: "); 
  for (int att=0; att < nAttributes; att++)
	printf("%c ", attIsContinuous[att]?'c':'d');
  printf("\n");
  printf("   hasmissingvalues: ");
  for (int att=0; att < nAttributes; att++)
	printf("%d ", hasMissingValues[att]);
  printf("\n");
  printf("attNValues (discrt): ");
  for (int att=0; att < nAttributes; att++)
	if (attIsContinuous[att]) printf("- ");
        else printf("%d ", attNValues[att]);
  printf("\n");
  printf("Target nValues: %d\n", targetNValues);

  // printf("Data: \n");
  for (int s=0; s < nSamples && s < nSamplesToPrint; s++)
	{
	for (int att=0; att < nAttributes; att++)
	   if (attIsContinuous[att])
		printf("%6.3f ", dataSample[s][att]);
	   else
		printf("%3d ", (int)dataSample[s][att]);
	printf("--> ");
	if (targetIsContinuous)
		printf("%5.3f", target[s]);
	else
		printf("%d", (int)target[s]);
	printf("\n");
	}

  if (nSamplesToPrint > 0 && nSamplesToPrint < nSamples)
	printf("(... %d samples omitted...)\n", nSamples-nSamplesToPrint);

  return;
  }


//-------------------------------------------------------

// reads a word, checks that it is the same as what we expected.
// fatal error if not.
void readWord(FILE *fp, char expectedWord[])
  {
  char buff[4096];

  safe_readStr(fp, buff);
  if (strcmp(buff, expectedWord))
	{
	fprintf(stderr, "ERROR: Expected %s, saw %s\n", expectedWord, buff);
	error("Error in readWord()");
	}

  return;
  }

// returns a double in general, or missingDataVal if "?"
double safe_readElement(FILE *fp, double missingDataVal)
  {
  char buff[4096];

  safe_readStr(fp, buff);
  if (!strcmp(buff, "?"))
	return missingDataVal;

  return safe_strtod(buff);
  }

// returns an integer from 0 to discretizationSize
int discretize(double x, double maxX, double minX, int discretizationSize)
  {
  assert(maxX >= x && x >= minX);
  if (maxX == minX)
	return 0;

  int d = (int)floor(discretizationSize*(x-minX)/(maxX-minX));
  if (d == discretizationSize)
	{
	assert(x==maxX);
	d = discretizationSize-1;
	}
  assert(d >= 0 && d < discretizationSize);

  return d;
  }

//-------------------------------------------------------

void readDataInfo(const char dataName[], int *nSamples, int *nAttributes, 
		Array<int> &attNValues, Array<int> &hasMissingValues,
		Array<int> &isDiscrete, 
		Array<double> &maxVals, Array<double> &minVals, int *targetIndex)
  {
  char infFilename[4096];
  char discreteContinuous[16384];

  strcpy(infFilename, DATADIR);		// ".../data/"
  strcat(infFilename, dataName); 	// ".../data/foo"
  strcat(infFilename, "/"); 	
  strcat(infFilename, dataName); 	
  strcat(infFilename, ".inf"); 		// ".../data/foo/foo.inf"

  FILE *fp = safe_fopen(infFilename, "rt");

  readWord(fp, "numberofsamples");
  *nSamples = safe_readInt(fp);

  readWord(fp, "numberofattributes");
  *nAttributes = safe_readInt(fp);

  readWord(fp, "targetattribute");
  *targetIndex = safe_readInt(fp);
  assert(*targetIndex >= 0 && *targetIndex < *nAttributes);

  attNValues.resize(*nAttributes);
  hasMissingValues.resize(*nAttributes);
  isDiscrete.resize(*nAttributes);
  maxVals.resize(*nAttributes); maxVals.setAllTo(0);
  minVals.resize(*nAttributes); minVals.setAllTo(0);

  readWord(fp, "discretecontinuous");
  assert(*nAttributes < 16384);
  safe_readStr(fp, discreteContinuous);
  assert(strlen(discreteContinuous) == *nAttributes);
  for (int att=0; att < *nAttributes; att++)
	{
	if (discreteContinuous[att] == 'd')
	    isDiscrete[att] = 1;
	else if (discreteContinuous[att] == 'c')
	    isDiscrete[att] = 0;
	else
	    error("Bad character in discrete/continuous info? (77251)");
	}

  for (int att=0; att < *nAttributes; att++)
	{
	readWord(fp, "attribute");
	int att0 = safe_readInt(fp);
	assert(att0 == att);
	readWord(fp, "nmissing");
	int nmissing = safe_readInt(fp);
	if (nmissing > 0)
		hasMissingValues[att] = 1;
	else 	
		hasMissingValues[att] = 0;

	if (isDiscrete[att])
	    {
	    char buff[4096];
	    readWord(fp, "nvalues");
	    int nValues = safe_readInt(fp);
	    attNValues[att] = nValues;
	    readWord(fp, "values");
	    for (int v=0; v < nValues; v++)
		    safe_readStr(fp, buff);	// skip the list of values
	    }
	else
	    {
	    attNValues[att] = INVALID_VALUE;
	    readWord(fp, "max");
	    maxVals[att] = safe_readDouble(fp);
	    readWord(fp, "min");
	    minVals[att] = safe_readDouble(fp);
	    }
	}

  fclose(fp);

  return;
  }


// returns 1 if this sample had 1 or more missing attributes.
int readSample(FILE *fp, Data *data, int nAttributes, int targetIndex, 
			Array<int> &hasMissingValues, int thisSampleIndex)
  {
  assert(hasMissingValues.length() == nAttributes);
  assert(data->dataSample.ncols() == nAttributes-1);

  const double MISSING_VALUE_TOKEN = data->missingValueToken();

  // dex indices into data->dataSample[thisSampleIndex]
  int dex=0;
  int thisHasMissingAtt = 0;
  for (int att=0; att < nAttributes; att++)
	{
	double val = safe_readElement(fp, MISSING_VALUE_TOKEN);
	if (val == MISSING_VALUE_TOKEN)
	  {
	  assert(hasMissingValues[att]);
	  thisHasMissingAtt = 1;
	  }

	if (att == targetIndex)
	    {
	    assert(val != MISSING_VALUE_TOKEN);	    // target cannot be missing-value
	    data->target[thisSampleIndex] = val;
	    }
	else	// not target value
	    {
	    data->dataSample[thisSampleIndex][dex] = val;
	    dex++;
	    }
	}
  assert(dex == nAttributes-1);

  return thisHasMissingAtt;
  }

void new_trainTestSplit(const Data &data,
			Data **trainData_p, Data **testData_p,
			int nTest)
  {
  Data *trainData, *testData;
  int nTrain;
  Array<int> assignToTest(data.nSamples);

  // select_frac_bits(assignToTest.allEle(), data.nSamples, testFract);
  // nTest = sum_arr(assignToTest.allEle(), data.nSamples);
  // nTrain = data.nSamples - nTest;

  // nTest = round(data.nSamples * testFract);
  // nTrain = data.nSamples - nTest;
  // select_nbits(assignToTest.allEle(), data.nSamples, nTest);

  assert(0 <= nTest && nTest <= data.nSamples);
  nTrain = data.nSamples - nTest;
  select_nbits(assignToTest.allEle(), data.nSamples, nTest);

  assert(data.nSamples >= 2);
  assert(nTrain >= 0 && nTest >= 0);
  if (nTrain == 0)
  	error("0 training samples in readDiscretized Data. Not good/not quite implemented.");
  // if (nTest == 0 && testFract != 0.0)
  //	error("0 test samples with non-zero test fraction. Not good....");

  trainData = new Data(nTrain, data);
  testData = new Data(nTest, data);

  int trainDex=0, testDex=0;
  for (int s=0; s < data.nSamples; s++)
    {
    if (assignToTest[s])
	{
	for (int att=0; att < data.nAttributes; att++)
	    testData->dataSample[testDex][att] = data.dataSample[s][att];
	testData->target[testDex] = data.target[s];
	testDex++;
	}
    else
	{
	for (int att=0; att < data.nAttributes; att++)
	    trainData->dataSample[trainDex][att] = data.dataSample[s][att];
	trainData->target[trainDex] = data.target[s];
	trainDex++;
	}
    }
  assert(trainDex == nTrain && testDex == nTest);

  trainData->shuffleSample();
  testData->shuffleSample();

  *trainData_p = trainData;
  *testData_p = testData;

  return;
  }

void dataSummary(const char dataName[], int *nSamples, int *nAttributes)
  {
  int targetIndex;
  Array<int> attNValues(1), hasMissingValues(1), isDiscrete(1);
  Array<double>	maxVals(1), minVals(1);

  readDataInfo(dataName, nSamples, nAttributes, attNValues, hasMissingValues,
		isDiscrete, maxVals, minVals, &targetIndex);

  *nAttributes -= 1;		// exclude target attribute

  return;
  }

// discardMissing means any data with 1 or more missing attributes will be
// discarded.  Note that this may make the dataset returned inconsistent
// with what dataSummary() returns.
Data* new_readData(const char dataName[], int discardMissing)
  {
  Data *data;
  int nSamples, nAttributes, targetIndex;
  Array<int> attNValues(1), hasMissingValues(1), isDiscrete(1);
  Array<double>	maxVals(1), minVals(1);

  readDataInfo(dataName, &nSamples, &nAttributes, attNValues, hasMissingValues,
		isDiscrete, maxVals, minVals, &targetIndex);

  // nAttributes counts target as an attribute as well.
  // for class Data, we'll need to subtract 1 back off for data->nAttributes;

  // allocate space for and set up internals of data
  data = new Data(nSamples, nAttributes-1);
			// nAttributes-1 because class Data does not count
			// target as an attribute, whereas our file format does.
  assert(data != NULL);

  int dex=0;
  for (int att=0; att < nAttributes; att++)
    {
    if (att == targetIndex)
	{
	data->targetIsContinuous = !isDiscrete[att];
	assert(!hasMissingValues[att]);		// missing values not allowed in target
	// printf("[%d]", hasMissingValues[att]);
	data->targetNValues = attNValues[att];
	}
    else
	{
	data->attIsContinuous[dex] = !isDiscrete[att];
	data->hasMissingValues[dex] = hasMissingValues[att];
	// printf("(%d)", hasMissingValues[att]);
	data->attNValues[dex] = attNValues[att];
	dex++;
	}
    }
  assert(dex == nAttributes-1);

  // start reading the data
  char datFilename[4096];
  strcpy(datFilename, DATADIR); strcat(datFilename, dataName); 	
  strcat(datFilename, "/"); strcat(datFilename, dataName); 	
  strcat(datFilename, ".dat"); 		// ".../data/foo/foo.dat"

  FILE *fp = safe_fopen(datFilename, "rt");

  int sampleDex=0, nHasMissingAtt=0;
  for (int ctr=0; ctr < nSamples; ctr++)
    {
    int hasMissingAtt = readSample(fp, data, nAttributes, targetIndex, hasMissingValues, sampleDex);
    if (discardMissing && hasMissingAtt)
      nHasMissingAtt++;
    else
      sampleDex++;
    }
  fclose(fp);

  // If some samples were discarded, then shrink the datastructures to reflect 
  // that we didn't actually get nSamples examples read in.
  if (!discardMissing)
    assert(nHasMissingAtt==0 && sampleDex==nSamples);
  if (discardMissing)
    {
    data->nSamples = nSamples - nHasMissingAtt;
    if (data->nSamples == 0)
      error("new_readData: After discarding missing, have no data left!");
    for (int dex=0; dex < data->nAttributes; dex++)
      data->hasMissingValues[dex] = 0;

    Matrix temp(data->nSamples, data->nAttributes);
    for (int i=0; i < data->nSamples; i++) for (int j=0; j < data->nAttributes; j++)
      temp[i][j] = data->dataSample[i][j];

    data->dataSample.resize(data->nSamples, data->nAttributes);
    data->dataSample = temp;
    }

  data->shuffleSample();

  assert(data->nAttributes == nAttributes-1);

  return data;
  }

Data* new_recodeMultivariateToBinary(const Data &oldData, int addIntercept)
  {
  // Calculate number of attributes in re-coded data
  int newNAttributes=0;
  for (int a=0; a < oldData.nAttributes; a++)
    if (oldData.attIsContinuous[a])
      newNAttributes++;
    else
      {
      assert(oldData.attNValues[a] >= 2);
      newNAttributes += oldData.attNValues[a]-1;
      }
  if (addIntercept)
    newNAttributes++;

  // Initialize summary datastructures for recoded data
  Data *newData = new Data(oldData.nSamples, newNAttributes);
  newData->targetIsContinuous = oldData.targetIsContinuous;
  newData->target = oldData.target;
  newData->targetNValues = oldData.targetNValues;

  int newADex=0;
  for (int a=0; a < oldData.nAttributes; a++)
    {
    if (oldData.attIsContinuous[a])
      {
      newData->attIsContinuous[newADex] = 1;
      newData->hasMissingValues[newADex] = oldData.hasMissingValues[a];
      assert(oldData.attNValues[a] == INVALID_VALUE);
      newData->attNValues[newADex] = INVALID_VALUE;
      newADex++;
      }
    else
      {
      for (int v=1; v < oldData.attNValues[a]; v++)
	{
	newData->attIsContinuous[newADex] = 0;
	newData->hasMissingValues[newADex] = oldData.hasMissingValues[a];
	newData->attNValues[newADex] = 2;
	newADex++;
	}
      }
    }
  if (addIntercept)
    {
    newData->attIsContinuous[newADex] = 1;
    newData->hasMissingValues[newADex] = 0;
    newData->attNValues[newADex] = INVALID_VALUE;
    newADex++;
    }
  assert(newADex == newNAttributes);

  // Finally, take care of the recoded data sample
  for (int ctr=0; ctr < oldData.nSamples; ctr++)
    {
    newADex=0;
    for (int a=0; a < oldData.nAttributes; a++)
      {
      if (oldData.attIsContinuous[a])
	{
	newData->dataSample[ctr][newADex] = oldData.dataSample[ctr][a];
	newADex++;
	}
      else
        {
        assert(isInt(oldData.dataSample[ctr][a]) 
	    || oldData.dataSample[ctr][a] == Data::missingValueToken());
        for (int v=1; v < oldData.attNValues[a]; v++)
  	  {
	  if (oldData.dataSample[ctr][a] == Data::missingValueToken())
	    newData->dataSample[ctr][newADex] = Data::missingValueToken();
	  else if ((int)oldData.dataSample[ctr][a] == v)
	    newData->dataSample[ctr][newADex] = 1;
	  else
	    newData->dataSample[ctr][newADex] = 0;
	  newADex++;
	  }
	}
      }
    if (addIntercept)
      {
      newData->dataSample[ctr][newADex] = 1;
      newADex++;
      }
    assert(newADex == newNAttributes);
    }

  return newData;
  }

Data* new_addIntercept(const Data &oldData)
  {
  // Initialize summary datastructures for new data 
  Data *newData = new Data(oldData.nSamples, oldData.nAttributes+1);
  newData->targetIsContinuous = oldData.targetIsContinuous;
  newData->target = oldData.target;
  newData->targetNValues = oldData.targetNValues;

  int newADex=0;
  for (int a=0; a < oldData.nAttributes; a++)
    {
    newData->attIsContinuous[a] = oldData.attIsContinuous[a];
    newData->hasMissingValues[a] = oldData.hasMissingValues[a];
    newData->attNValues[a] = oldData.attNValues[a];
    }
  newData->attIsContinuous[oldData.nAttributes] = 1;
  newData->hasMissingValues[oldData.nAttributes] = 0;
  newData->attNValues[oldData.nAttributes] = INVALID_VALUE;

  // Finally, take care of the recoded data sample
  for (int ctr=0; ctr < oldData.nSamples; ctr++)
    {
    for (int a=0; a < oldData.nAttributes; a++)
      newData->dataSample[ctr][a] = oldData.dataSample[ctr][a];
    newData->dataSample[ctr][oldData.nAttributes] = 1;
    }

  return newData;
  }

