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
  double lower_bound = 0;
  double upper_bound = 100;
  std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
  std::default_random_engine re;

  std::vector<double> vals;
  std::vector<int> offsets;
  std::vector<int> cols;

  int numTotal = 0;
  for (int r=0; r < nRows; ++r){
    offsets.push_back(numTotal);
    int numInRow = 0;
    //printf("[");
    for (int c=0; c < nCols; ++c){
      double test = unif(re);
      if (test < sparsePercent){
        double val = unif(re);
        vals.push_back(val);
        ++numInRow;
        cols.push_back(c);
        //printf(" %2d:%8.4f", c, val);
      }
    }
    //printf(" ]\n");
    numTotal += numInRow;
  }
  offsets.push_back(numTotal);
  std::cout << "No. nonzeros: " << numTotal << std::endl;

  CSRMatrix<double> mat(std::move(vals), std::move(offsets), std::move(cols), nCols);

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
