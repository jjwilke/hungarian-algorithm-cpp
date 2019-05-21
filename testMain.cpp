#include <iostream>
#include "Hungarian.h"


int main(void)
{
  //M[i + nRows * j] = W(i,j)
  //each column is the weights of a given vertex to every other vertex
  vector<double> costMatrix = { 10, 19, 8, 15, 
                                10, 18, 7, 17,
                                13, 16, 9, 14, 
                                12, 19, 8, 18 };

	HungarianAlgorithm HungAlgo;
	vector<int> assignment;

  int nRows = 4;
  int nCols = 4;
  double cost = HungAlgo.Solve(costMatrix, assignment, nRows, nCols);

  double costCheck = 0;
	for (unsigned int x = 0; x < nRows; x++){
		std::cout << x << "," << assignment[x] << "\t";
    int edgeIndex = x + nRows*assignment[x];
    costCheck += costMatrix[edgeIndex];
  }

	std::cout << "\ncost:  " << cost << std::endl;
  std::cout <<  "check:  " << costCheck << std::endl;

	return 0;
}
