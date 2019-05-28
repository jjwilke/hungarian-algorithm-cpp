///////////////////////////////////////////////////////////////////////////////
// Hungarian.cpp: Implementation file for Class HungarianAlgorithm.
// 
// This is a C++ wrapper with slight modification of a hungarian algorithm implementation by Markus Buehren.
// The original implementation is a few mex-functions for use in MATLAB, found here:
// http://www.mathworks.com/matlabcentral/fileexchange/6543-functions-for-the-rectangular-assignment-problem
// 
// Both this code and the orignal code are published under the BSD license.
// by Cong Ma, 2016
// 

#include <stdlib.h>
#include <cfloat> // for DBL_MAX
#include <cmath>  // for fabs()
#include "Hungarian.h"

#define IS_ZERO(x) (fabs(x) < DBL_EPSILON)

//#define debugLine(str, ...) printf(str "\n", __VA_ARGS__)
//#define debug(str, ...) printf(str, __VA_ARGS__)
//#define debugNL() printf("\n"); fflush(stdout)

#define debugLine(str, ...)
#define debug(str, ...)
#define debugNL()

HungarianAlgorithm::HungarianAlgorithm(const CSRMatrix<double>& matrix) :
  coveredColumns_(matrix.nCols(), -1),
  coveredRows_(matrix.nRows(), -1),
  starredCols_(matrix.nCols(), -1),
  starredRows_(matrix.nRows(), -1)
{
}

//********************************************************//
// A single function wrapper for solving assignment problem.
//********************************************************//
// Mind the index is "i + nRows * j".
// Here the cost matrix of size MxN is defined as a double precision array of N*M elements. 
double HungarianAlgorithm::solve(const CSRMatrix<double>& matrix)
{
  //we are programmed to find the min...
  //'invert' the elements to make this a max problm
  double max = matrix.max();
  auto inverted = matrix.mutate([=](double inp) -> double { return max + 1 - inp; } );

  std::fill(starredCols_.begin(), starredCols_.end(), -1);
  std::fill(starredRows_.begin(), starredRows_.end(), -1);
  std::fill(coveredColumns_.begin(), coveredColumns_.end(), false);
  std::fill(coveredRows_.begin(), coveredRows_.end(), false);

	// call solving function
  step2(inverted);
  return computeCost(matrix);
}


//********************************************************//
// Solve optimal solution for assignment problem using Munkres algorithm, also known as Hungarian Algorithm.
//********************************************************//
void HungarianAlgorithm::findOptimal(CSRMatrix<double>& matrix)
{
  /* move to step 2b */
  step2(matrix);
}



/********************************************************/
double HungarianAlgorithm::computeCost(const CSRMatrix<double>& matrix)
{
  double cost = 0;
  int nOfRows = matrix.nRows();
  for (int row = 0; row<nOfRows; row++){
    //assignment stores RELATIVE col indices
    int col = starredRows_[row];
    if (col >= 0){
      int rowOffset = matrix.rowOffset(row);
      cost += matrix.valueAt(rowOffset + col);
    }
	}
  return cost;
}

/********************************************************/
void HungarianAlgorithm::step2(CSRMatrix<double>& matrix)
{
  /* preliminary steps */
  int nOfRows = matrix.nRows();
  double* rowValues;
  int* colInds;
  int nCols;
  std::vector<double> colMins(matrix.nCols(), std::numeric_limits<double>::max());
  for (int row = 0; row<nOfRows; row++){
    std::tie(nCols, rowValues, colInds) = matrix.getRow(row);
    double rowMin = std::numeric_limits<double>::max();
    //find the minimum
    for (int col=0; col < nCols; ++col){
      double val = rowValues[col];
      rowMin = std::min(rowMin, val);
    }
    //substract the minimum in a normalization pass
    for (int col=0; col < nCols; ++col){
      rowValues[col] -= rowMin;
      int absCol = colInds[col];
      colMins[absCol] = std::min(rowValues[col], colMins[absCol]);
    }
  }

  //we need another pass to substract column minima and find zeros
  for (int row=0; row < nOfRows; ++row){
    std::tie(nCols, rowValues, colInds) = matrix.getRow(row);
    for (int col=0; col < nCols; ++col){
      int absCol = colInds[col];
      rowValues[col] -= colMins[absCol];
      debug(" %d:%8.4f", absCol, rowValues[col]);
    }
    debugNL();
    //find the first zero, set booleans and break
  }


	/* move to step 3 */
  step3(matrix);
}


/********************************************************/
void HungarianAlgorithm::step3(CSRMatrix<double>& matrix)
{
  int numSteps = 0;
  while (true){
    //this is the "drawing" portion where we have to mark all rows and columns with zeros
    std::fill(coveredRows_.begin(), coveredRows_.end(), false);
    std::fill(coveredColumns_.begin(), coveredColumns_.end(), false);
    std::fill(starredCols_.begin(), starredCols_.end(), -1);
    int numOfRows = matrix.nRows();
    double* rowValues;
    int* colInds;
    int nCols;
    int numStarredRows = 0;
    int numCoveredRows = 0;
    int numCoveredCols = 0;

    for (int row=0; row < numOfRows; ++row){
      std::tie(nCols, rowValues, colInds) = matrix.getRow(row);
      starredRows_[row] = -1;
      for (int col=0; col < nCols; ++col){
        int absCol = colInds[col];
        if (IS_ZERO(rowValues[col]) && starredCols_[absCol] == -1){
          debugLine("Starring and covering row=%d col=%d:%d", row, absCol, col);
          starredCols_[absCol]  = row;
          starredRows_[row] = col;
          coveredRows_[row] = true;
          ++numStarredRows;
          ++numCoveredRows;
          break;
          //otherwise keep going until I find another zero
        }
      }
    }

    if (numStarredRows == numOfRows){
      return; //we got em all!
    } else {
      debugLine("Continuing to later steps with %d covered rows", numStarredRows);
    }

    int nextRow = 0;
    int row = 0;
    std::vector<int> uncoveredRows;
    while (nextRow < numOfRows || !uncoveredRows.empty()){
      //we are either looping through uncovered rows
      //or bouncing around to columns

      if (uncoveredRows.empty()){
        row = nextRow++;
      } else {
        row = uncoveredRows.back();
        uncoveredRows.pop_back();
      }
      if (!coveredRows_[row]){
        debugLine("Visting row %d", row);
        std::tie(nCols, rowValues, colInds) = matrix.getRow(row);
        for (int col=0; col < nCols; ++col){
          int absCol = colInds[col];
          if (row == 0){
            debugLine("Row[%d:%d] = %8.4f ? %s", col, absCol, rowValues[col],
                      (coveredColumns_[absCol] ? "true" : "false"));
          }
          if (IS_ZERO(rowValues[col]) && !coveredColumns_[absCol]){
            if (starredCols_[absCol] != -1){ //cover the zeroes
              //cover the column, uncover the row previously covered
              int uncoveredRow = starredCols_[absCol];
              debugLine("Covering column %d, uncovering row %d", absCol, uncoveredRow);
              coveredRows_[uncoveredRow] = false;
              --numCoveredRows;
              uncoveredRows.push_back(uncoveredRow);
            } else {
              debugLine("Just covering column %d", absCol);
            }
            coveredColumns_[absCol] = true;
            ++numCoveredCols;
          }
        }
      }
    }

    if ((numCoveredCols + numCoveredRows) == numOfRows){
      std::cerr << "I have a unique solution!" << std::endl;
      abort();
    }


    /**
      while (!rowsToVisit.empty()){
        int nextRow = rowsToVisit.front();
        rowsToVisit.pop_front();
        visitRowStep3(matrix, nextRow, rowsToVisit, numCoveredRows, numCoveredCols);
      }
      if (coveredRows_[row] == -1){
        debugLine("Row %d is unassigned", row);
        visitRowStep3(matrix, row, rowsToVisit, numCoveredRows, numCoveredCols);
      }
    }
    //clean up deque at the end
    while (!rowsToVisit.empty()){
      int nextRow = rowsToVisit.front();
      rowsToVisit.pop_front();
      visitRowStep3(matrix, nextRow, rowsToVisit, numCoveredRows, numCoveredCols);
    }
    */

    step4(matrix);
    ++numSteps;
  }
}

void HungarianAlgorithm::findKnownOptimal(CSRMatrix<double>& matrix)
{
  /**
  int numOfRows = matrix.nRows();
  double* rowValues;
  int* colInds;
  int nCols;
  bool newZeros = true;
  int numCoveredRows = 0;
  while (newZeros && numCoveredRows < numOfRows){
    newZeros = false;
    numCoveredRows = 0;
    for (int row=0; row < numOfRows; ++row){
      std::tie(nCols, rowValues, colInds) = matrix.getRow(row);
      //try to move this to another column, if possible
      for (int col=0; col < nCols; ++col){
        int absCol = colInds[col];
        if (IS_ZERO(rowValues[col]) && coveredColumns_[absCol] == -1){
          newZeros = true;
          if (coveredRows_[row] != -1){
            int oldCol = colInds[coveredRows_[col]];
            coveredColumns_[oldCol] = -1;
            debugLine("Changing row=%d to col=%d:%d", row, absCol, col);
          } else {
            debugLine("First cover row=%d to col=%d:%d", row, absCol, col);
          }
          ++numCoveredRows;
          coveredRows_[row] = col;
        }
      }
    }
  }
  */
}

void HungarianAlgorithm::visitRowStep3(CSRMatrix<double>& matrix, int row, std::deque<int>& rowsToVisit,
                                       int& numCoveredRows, int& numCoveredCols)
{
  /**
  //I don't have an assignment
  markedRows_[row] = true;
  --numCoveredRows;
  double* rowValues;
  int* colInds;
  int nCols;
  std::tie(nCols, rowValues, colInds) = matrix.getRow(row);
  for (int col=0; col < nCols; ++col){
    if (IS_ZERO(rowValues[col])){
      int absCol = colInds[col];
      if (!markedCols_[absCol]){
        debugLine("\tCol %d is unmarked", absCol);
        markedCols_[absCol] = true;
        ++numCoveredCols;
        if (coveredColumns_[absCol] != -1){ //this column is assigned to someone
          int rowToMark = coveredColumns_[absCol]; //mark them
          if (coveredRows_[rowToMark] != -1 && !markedRows_[rowToMark]){
            debugLine("\t\tMarking assigned row %d", rowToMark);
            markedRows_[rowToMark] = true;
            rowsToVisit.push_back(rowToMark);
          }
        }
      }
    }
  }
  */
}

/********************************************************/
void HungarianAlgorithm::step4(CSRMatrix<double>& matrix)
{
  //find the minimum in the unmarked portion
  double h = std::numeric_limits<double>::max();
  int nOfRows = matrix.nRows();
  double* rowValues;
  int* colInds;
  int nCols;

  for (int row=0; row < nOfRows; ++row){
    if (coveredRows_[row]){
      debugLine("Row %d covered", row);
    }
  }

  for (int col=0; col < matrix.nCols(); ++col){
    if (coveredColumns_[col]){
      debugLine("Cols %d covered", col);
    }
  }

  //the "covering" of zeros is all marked columns and unmarked rows
  for (int row = 0; row<nOfRows; row++){
    if (!coveredRows_[row]){
      std::tie(nCols, rowValues, colInds) = matrix.getRow(row);
      for (int col = 0; col< nCols; col++){
        int absCol = colInds[col];
        if (!coveredColumns_[absCol]){
          h = std::min(h, rowValues[col]);
        }
      }
    }
  }
  if (h == std::numeric_limits<double>::max()){
    std::cerr << "No minimum found - nothing covered!" << std::endl;
    abort();
  }
  if (h == 0){
    std::cerr << "Zero minimum found - left zeroes uncovered!" << std::endl;
    abort();
  }
  debugLine("Substracting minimum %10.4f", h);

  /* add h to each unmarked row, subtract h from each uncovered column
   * the net effect is to substract from all uncovered, add to all covered twice
  */
  for (int row = 0; row<nOfRows; row++){
    double rowAdd = coveredRows_[row] ? h : 0;
    std::tie(nCols, rowValues, colInds) = matrix.getRow(row);
    for (int col = 0; col<nCols; col++){
      int absCol = colInds[col];
      double colAdd = coveredColumns_[absCol] ? 0 : -h;
      rowValues[col] += rowAdd + colAdd;
      debug(" %d:%d %8.4f", col, absCol, rowValues[col]);
    }
    debugNL();
  }
}


