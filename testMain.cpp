#include <iostream>
#include "Hungarian.h"
#include <random>


int main(int argc, char** argv)
{
  //M[i + nRows * j] = W(i,j)
  //each column is the weights of a given vertex to every other vertex
  //vector<double> costMatrix = { 10, 19, 8, 15, 
  //                            10, 18, 7, 17,
  //                            13, 16, 9, 14, 
  //                            12, 19, 8, 18 };

  int nRows = atoi(argv[1]);
  int sparsePercent = atoi(argv[2]);
  int nCols = nRows;
  std::vector<double> costMatrix(nRows*nCols);
  double lower_bound = 0;
  double upper_bound = 100;
  std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
  std::default_random_engine re;
  for (int r=0; r < nRows; ++r){
    printf("[");
    for (int c=0; c < nCols; ++c){
      int edgeIndex = r + c*nRows;
      double test = unif(re);
      if (test < sparsePercent){
        costMatrix[edgeIndex] = 0;
      } else {
        costMatrix[edgeIndex] = unif(re);
        //costMatrix[edgeIndex] = (r*nCols + c) % 10;
      }
      printf(" %8.4f", costMatrix[edgeIndex]);
    }
    printf(" ]\n");
  }

  auto mat = CSRMatrix<double>::createDense(std::move(costMatrix), nRows, nCols);

  HungarianAlgorithm HungAlgo(mat);
  double cost = HungAlgo.solve(mat);

  auto& assignment = HungAlgo.assignment();

  double costCheck = 0;
  for (int row = 0; row < nRows; row++){
    int rowOffset = mat.rowOffset(row);
    int col = assignment[row];
    int absCol = mat.colAt(rowOffset + col);
    std::cout << row << "," << absCol << "\t";
    costCheck += mat.valueAt(rowOffset + col);
  }

	std::cout << "\ncost:  " << cost << std::endl;
  std::cout <<  "check:  " << costCheck << std::endl;

	return 0;
}
