#include <stdio.h>
#include <values.h>

#include "assocarr.h"
#include "dym_arr.h"
#include "str.h"
#include "matrix.h"
#include "misc.h"

struct DiscreteVariable
  {
  private:
	String name;
	int nvalues;
	String *values;

  public:
	DiscreteVariable(const String &name_, int nvalues_, const String *values_);
	~DiscreteVariable();
	void setName(const String &newName);
	void setValues(int nvalues_, const String *values_);

	const String &nameVal(void) { return name; }
	int nvaluesVal(void) { return nvalues; }
	const String *valuesVal(void) { return values; }
  };

DiscreteVariable::DiscreteVariable(const String &name_, int nvalues_, 
				const String *values_)
	: name(name_), nvalues(nvalues_)
  {
  values = new String[nvalues];
  for (int ctr=0; ctr < nvalues; ctr++)
  	values[ctr] = values_[ctr];
  }

DiscreteVariable::~DiscreteVariable()
  {
  delete[] values;
  }

void DiscreteVariable::setName(const String &newName)
  { 
  name=newName; 
  }

void DiscreteVariable::setValues(int nvalues_, const String *values_)
  {
  delete[] values;

  nvalues = nvalues_;
  values = new String[nvalues];
  for (int ctr=0; ctr < nvalues; ctr++)
  	values[ctr] = values_[ctr];
  }

//------------------------------------------------------------

class JunctionTreeClique
  {
  public:
	int nvariables;
	DiscreteVariable *variables;
  };

class JunctionTreeEdge
  {
  int node1, node2;		// might be interpreted as a directed arc
				// from node1 to node2, or not. Someone else keeps
				// track;

  public:
	JunctionTreeEdge(int node1_, int node2_) : node1(node1_), node2(node2_) {}
  };

class JunctionTree
  {
  private:
	DiscreteVariable* variables;
	AssocArray<int,int> *parent;	// in the directed-graph interpretation,
					// parent[i][j] is 1 if j is a parent of i
	int nvariables;
	int graphIsDirected;


	void checkAdjacencyList(void);
	int connected(int i, int j) { return parent[i].query_val(j)
					    || parent[j].query_val(i); }
	int nunjoinedNeighborsOfNode(int i);
	void findCliques(const int decimationOrder[]);
	void dpinToUpin(void);
	int isSubcliqueOf(int *newClique, int newCliqueSize, HashedSet<int> &clique);
	void triangulateGraph(int decimationOrder[]);

  public:
	void addDirectedArc(const String &fromNode, const String &toNode);
	void dpinToJunctionTree(void);

	void print(void);

	JunctionTree(int nvariables_);
	~JunctionTree();
  };

JunctionTree::JunctionTree(int nvariables_) 
		: nvariables(nvariables_), graphIsDirected(1)
  {
  String defaultValues[2] = {"False", "True"};

  variables = new DiscreteVariable[nvariables]("", 2, defaultValues);
  for (int ctr=0; ctr < nvariables; ctr++)
	variables[ctr].setName("Variable" + as_string(ctr));

  parent = new AssocArray<int,int>[nvariables](0);

  return;
  }

JunctionTree::~JunctionTree()
  {
  delete[] variables;
  delete[] parent;

  return;
  }

void JunctionTree::dpinToUpin(void)
  {
  Dymarr<int> edgesToAdd(-1);
  int edgesToAddCtr = 0;
  int *myParents;

  if (!graphIsDirected)
	error("JunctionTree::dpinToUpin called on undirected graph");

  // Step 1: moralize parents

  // consecutive (odd,even-indexed) pairs on edgesToAdd will represent
  // links we have to add to the graph
  for (int ctr=0; ctr < nvariables; ctr++)
    {
    int nele = parent[ctr].nele_val();
    myParents = parent[ctr].new_all_keys();

    for (int dex1=0; dex1 < nele; dex1++)
      for (int dex2=dex1+1; dex2 < nele; dex2++)
	{
	int var1 = myParents[dex1];
	int var2 = myParents[dex2];
	if (!connected(var1, var2))
		{
		edgesToAdd[edgesToAddCtr++] = var1;
		edgesToAdd[edgesToAddCtr++] = var2;
		}
	}
    assert((myParents != NULL) ^ (nele == 0));
    if (myParents != NULL)
	delete[] myParents;
    }

  // actually add the edges to the graph
  assert(is_even(edgesToAddCtr));
  for (int ctr=0; ctr < edgesToAddCtr; ctr+=2)
    {
    assert(parent[edgesToAdd[ctr]][edgesToAdd[ctr+1]] == 0);
    parent[edgesToAdd[ctr]][edgesToAdd[ctr+1]] = 1;
    }

  // Step 2: do all the adjacency lists so that a->b implies b->a

  for (int var1=0; var1 < nvariables; var1++)
    {
    int nele = parent[var1].nele_val();
    myParents = parent[var1].new_all_keys();
    for (int ctr=0; ctr < nele; ctr++)
	{
	int var2 = myParents[ctr];
	assert(parent[var1][var2] == 1);
	parent[var2][var1] = 1;
	}
    assert((myParents != NULL) ^ (nele == 0));
    if (myParents != NULL)
	delete[] myParents;
    }

  graphIsDirected = 0;

  checkAdjacencyList();

  return;
  }

// very inefficient printing algorithm.
void JunctionTree::print(void)
  {
  printf("Variables: ");
  for (int ctr=0; ctr < nvariables; ctr++)
	{
	printf("%s ", variables[ctr].nameVal().as_char());
	if ((ctr+1) % 4 == 0 && ctr != nvariables-1)
		printf("\n           ");
	}
  printf("\n\n");

  printf("Edges:\n");
  for (int var1=0; var1 < nvariables; var1++)
    for (int var2=0; var2 < nvariables; var2++)
	{
	if (var1 == var2)	
		continue;
	if (parent[var2].query_val(var1))
	   {
	   if (!graphIsDirected && var1 < var2)
		printf((variables[var1].nameVal() + " - " + variables[var2].nameVal() + "\n").as_char());
	   else if (graphIsDirected)
		printf((variables[var1].nameVal() + " -> " + variables[var2].nameVal() + "\n").as_char());
	   }
	}

  return;
  }

int JunctionTree::nunjoinedNeighborsOfNode(int i)
  {
  int nunjoinedNeighbors=0;
  int *neighbors, nele;

  assert(!graphIsDirected);

  // get the list of our neighbors
  nele = parent[i].nele_val();
  neighbors = parent[i].new_all_keys();

  // check each pair for connectedness
  for (int ctr1=0; ctr1 < nele; ctr1++)
    for (int ctr2=ctr1+1; ctr2 < nele; ctr2++)
	{
	if (!connected(neighbors[ctr1], neighbors[ctr2]))
		{
		nunjoinedNeighbors++;
		// printf("%d,%d unconnected for %d\n",
		//	neighbors[ctr1], neighbors[ctr2], i);
		}
	}

  if (neighbors != NULL)
	delete[] neighbors;

  return nunjoinedNeighbors;
  }

void JunctionTree::checkAdjacencyList(void)
  {
  for (int ctr=0; ctr < nvariables; ctr++)
	{
	int nele = parent[ctr].nele_val();
	int *neighbors = parent[ctr].new_all_keys();

	for (int dex=0; dex < nele; dex++)
	    {
	    assert(parent[ctr][neighbors[dex]] == 1);
	    if (!graphIsDirected)
		assert(parent[neighbors[dex]][ctr] == 1);
	    }

	if (neighbors != NULL)
		delete[] neighbors;
	}

  return;
  }

// determines if the elements of newClique is a subset of clique
int JunctionTree::isSubcliqueOf(int *newClique, int newCliqueSize, HashedSet<int> &clique)
  {
  for (int ctr=0; ctr < newCliqueSize; ctr++)
    if (!clique.containsElement(newClique[ctr]))
	return 0;
  return 1;
  }

// Precond: graph is triangulated, and decimationOrder[] gives
// a correct decimation order.
void JunctionTree::findCliques(const int decimationOrder[])
  {
  int *possibleNewClique = new int[nvariables];
  HashedSet<int> *cliques = new HashedSet<int>[nvariables];
  int ncliques = 0;
  int *decimated = new int[nvariables];
  zero_array(decimated, nvariables);

  for (int dex=0; dex < nvariables; dex++)
    {
    int nextToDecimate = decimationOrder[dex];

    // the to-be-decimated node and its neighbors form a potential clique
    int possibleNewCliqueSize = 0;
    int *neighbors = parent[nextToDecimate].new_all_keys();
    for (int ctr=0; ctr < parent[nextToDecimate].nele_val(); ctr++)
	{
	if (!decimated[neighbors[ctr]])
	    possibleNewClique[possibleNewCliqueSize++] = neighbors[ctr];
	assert(parent[nextToDecimate][neighbors[ctr]] == 1);
	}
    possibleNewClique[possibleNewCliqueSize++] = nextToDecimate;

    int isSubclique = 0;
    for (int ctr=0; ctr < ncliques && !isSubclique; ctr++)
	if (isSubcliqueOf(possibleNewClique, possibleNewCliqueSize, cliques[ctr]))
		{
		isSubclique = 1;
		printf("While decimating %d, found subclique of number %d\n",
			nextToDecimate, ctr);
		break;
		}
    if (!isSubclique)
	{
	printf("New clique while decimating %d: ", nextToDecimate);
	// We've found a new clique.
	for (int i=0; i < possibleNewCliqueSize; i++)
	    {
	    cliques[ncliques].insertElement(possibleNewClique[i]);
	    printf("%d ", possibleNewClique[i]);
	    }
	printf("\n");
	ncliques++;
	}

    decimated[nextToDecimate] = 1;

    if (neighbors != NULL)
	delete[] neighbors;
    }

  // HERE: after finding the cliques, we should probably build a junction
  // tree out of them.

  printf("Cliques: \n");
  for (int ctr=0; ctr < ncliques; ctr++)
    {
    printf("    ");
    int *cliqueElements = cliques[ctr].new_allElements();
    assert(cliqueElements != NULL && cliques[ctr].nele_val() > 0);
    for (int dex=0; dex < cliques[ctr].nele_val(); dex++)
	printf((variables[cliqueElements[dex]].nameVal()+" ").as_char());
    printf("\n");
    delete[] cliqueElements;
    }

  delete[] decimated;
  delete[] cliques;
  delete[] possibleNewClique;

  return;
  }

void JunctionTree::dpinToJunctionTree(void)
  {
  int *decimationOrder = new int[nvariables];

  // Moralize and convert to undirected
  dpinToUpin();

  // Triangulate
  triangulateGraph(decimationOrder);

  // Find cliques
  findCliques(decimationOrder);

  printf("------------------------------\n");
  printf("Triangulated: \n");
  print();

  delete[] decimationOrder;

  return;
  }

// precond: Graph has been made into a UPIN already
void JunctionTree::triangulateGraph(int decimationOrder[])
  {
  if (graphIsDirected)
	error("JunctionTree::triangulateGraph called on directed graph");

  checkAdjacencyList();

  int *nunjoinedNeighbors = new int[nvariables];
  int *decimated = new int[nvariables];
  for (int ctr=0; ctr < nvariables; ctr++)
	nunjoinedNeighbors[ctr] = nunjoinedNeighborsOfNode(ctr);
  zero_array(decimated, nvariables);

  for (int dex=0; dex < nvariables; dex++)
    {
    int nextToDecimate = argmin(nunjoinedNeighbors, nvariables);
    int edgesAdded = 0;

    decimationOrder[dex] = nextToDecimate;
    decimated[nextToDecimate] = 1;

    int nele = parent[nextToDecimate].nele_val();
    int *neighbors = parent[nextToDecimate].new_all_keys();

    printf("%d %d %d %d %d %d %d %d (%d)\n",
	nunjoinedNeighbors[0], nunjoinedNeighbors[1], 
	nunjoinedNeighbors[2], nunjoinedNeighbors[3], 
	nunjoinedNeighbors[4], nunjoinedNeighbors[5], 
	nunjoinedNeighbors[6], nunjoinedNeighbors[7], 
	nextToDecimate);

    // connect all unconnected neighbors 
    for (int ctr1=0; ctr1 < nele; ctr1++)
      {
      assert(neighbors[ctr1] != nextToDecimate);
      if (decimated[neighbors[ctr1]])
	  continue;
      for (int ctr2=ctr1+1; ctr2 < nele; ctr2++)
	{
        assert(neighbors[ctr2] != nextToDecimate);
	if (decimated[neighbors[ctr2]])
		continue;
	if (!connected(neighbors[ctr1], neighbors[ctr2]))
	    {
	    int old1nele = parent[neighbors[ctr1]].nele_val();
	    int *old1Neighbors = parent[neighbors[ctr1]].new_all_keys();
	    int old2nele = parent[neighbors[ctr2]].nele_val();
	    int *old2Neighbors = parent[neighbors[ctr2]].new_all_keys();
	    for (int i=0; i < old1nele; i++)
		{
		if (!decimated[old1Neighbors[i]]
		   	&& !connected(old1Neighbors[i], neighbors[ctr2]))
		    {
		    nunjoinedNeighbors[neighbors[ctr1]]++;
		    printf("%d, %d new unjoined for %d \n", 
			old1Neighbors[i], neighbors[ctr2], neighbors[ctr1]);
		    }
		if (!decimated[old1Neighbors[i]]
			&& parent[neighbors[ctr2]].query_val(old1Neighbors[i]))
		    {
		    // common neighbor
		    printf("Adding %d,%d deducts 1 for %d\n",
			neighbors[ctr1], neighbors[ctr2], old1Neighbors[i]);
		    nunjoinedNeighbors[old1Neighbors[i]]--;
		    }
		}
	    for (int i=0; i < old2nele; i++)
		if (!decimated[old2Neighbors[i]]
			&& !connected(old2Neighbors[i], neighbors[ctr1]))
		    {
		    nunjoinedNeighbors[neighbors[ctr2]]++;
		    printf("%d, %d new unjoined for %d\n", 
			old2Neighbors[i], neighbors[ctr1], neighbors[ctr2]);
		    }
	    if (old1Neighbors != NULL)
		delete[] old1Neighbors;
	    if (old2Neighbors != NULL)
		delete[] old2Neighbors;

	    parent[neighbors[ctr1]][neighbors[ctr2]] = 1;
	    parent[neighbors[ctr2]][neighbors[ctr1]] = 1;
	    printf("Adding edge: %d,%d\n", 
			neighbors[ctr1], neighbors[ctr2]);
	    edgesAdded++;
	    }

	}
      }

    assert(edgesAdded == nunjoinedNeighbors[nextToDecimate]);
    nunjoinedNeighbors[nextToDecimate] = MAXINT;

    // update my neighbor's nunjoinedNeighbors count 
    for (int ctr1=0; ctr1 < nele; ctr1++)
	{
	int myNeighbor = neighbors[ctr1];
	if (decimated[myNeighbor])	// note excludes neighbor[ctr1] as well
		continue;
	int *neighborsNeighbors = parent[myNeighbor].new_all_keys();
	int nele2 = parent[myNeighbor].nele_val();
	for (int ctr2=0; ctr2 < nele2; ctr2++)
	    {
	    if (decimated[neighborsNeighbors[ctr2]])
		continue;
	    assert(parent[myNeighbor][neighborsNeighbors[ctr2]]==1);
	    if (!connected(nextToDecimate, neighborsNeighbors[ctr2]))
		{
		printf("%d, %d no longer unjoined for %d\n",
			nextToDecimate, neighborsNeighbors[ctr2],
			myNeighbor);
		nunjoinedNeighbors[myNeighbor]--;
		assert(nunjoinedNeighbors[myNeighbor]>=0);
		}
	    }
	if (neighborsNeighbors != NULL)
		delete[] neighborsNeighbors;
	}

    if (neighbors != NULL)
	delete[] neighbors;
    }

  delete[] nunjoinedNeighbors;
  delete[] decimated;

  return;
  }

void JunctionTree::addDirectedArc(const String &fromNode, const String &toNode)
  {
  int fromDex= -1, 
	toDex= -1;

  if (!graphIsDirected)
	error("JunctionTree::addDirectedArc called on undirected graph: not good.");

  for (int ctr=0; ctr < nvariables; ctr++)
    if (fromNode == variables[ctr].nameVal())
	{ fromDex = ctr; break; }
  if (fromDex == -1)
	error(("JunctionTree::addDirectedArc: " + fromNode + " not found.").as_char());
		
  for (int ctr=0; ctr < nvariables; ctr++)
    if (toNode == variables[ctr].nameVal())
	{ toDex = ctr; break; }
  if (toDex == -1)
	error(("JunctionTree::addDirectedArc: " + toNode + " not found.").as_char());

  parent[toDex][fromDex] = 1;

  return;
  };

//-------------------------------------------------------

int main(void)
  {
  JunctionTree jt(8);

  jt.addDirectedArc("Variable0", "Variable1");
  jt.addDirectedArc("Variable0", "Variable2");
  jt.addDirectedArc("Variable1", "Variable3");
  jt.addDirectedArc("Variable2", "Variable4");
  jt.addDirectedArc("Variable3", "Variable5");
  jt.addDirectedArc("Variable4", "Variable5");

//  jt.addDirectedArc("Variable6", "Variable4");
//  jt.addDirectedArc("Variable6", "Variable5");
//  jt.addDirectedArc("Variable6", "Variable7");

  jt.addDirectedArc("Variable4", "Variable6");
  jt.addDirectedArc("Variable4", "Variable7");
  jt.addDirectedArc("Variable5", "Variable6");
  jt.addDirectedArc("Variable7", "Variable6");

  printf("------------------------------\n");
  printf("Original graph: \n");
  jt.print();

  jt.dpinToJunctionTree();

  };

