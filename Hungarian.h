///////////////////////////////////////////////////////////////////////////////
// Hungarian.h: Header file for Class HungarianAlgorithm.
// 
// This is a C++ wrapper with slight modification of a hungarian algorithm implementation by Markus Buehren.
// The original implementation is a few mex-functions for use in MATLAB, found here:
// http://www.mathworks.com/matlabcentral/fileexchange/6543-functions-for-the-rectangular-assignment-problem
// 
// Both this code and the orignal code are published under the BSD license.
// by Cong Ma, 2016
// 

#ifndef CSR_HUNGARIAN_H
#define CSR_HUNGARIAN_H

#include <iostream>
#include <vector>
#include <tuple>
#include <deque>

using namespace std;

template <class T>
class CSRMatrix
{
 public:
  CSRMatrix(std::vector<T>&& vals, std::vector<int>&& rowOffsets, std::vector<int>&& columns, int ncols)
    : values_(std::move(vals)), rowOffsets_(std::move(rowOffsets)), columns_(std::move(columns)),
      nCols_(ncols)
  {
  }

  int nRows() const {
    return rowOffsets_.size() - 1;
  }

  int nCols() const {
    return nCols_;
  }

  T max() const {
    T myMax = std::numeric_limits<T>::min();
    for (T val : values_){
      myMax = std::max(val, myMax);
    }
    return myMax;
  }

  static CSRMatrix createDense(std::vector<T>&& vals, int nRows, int nCols){
    std::vector<int> rowOffsets(nRows+1);
    std::vector<int> colInds(nRows*nCols);
    for (int i=0; i < nRows; ++i){
      rowOffsets[i] = i*nCols;
      for (int j=0; j < nCols; ++j){
        colInds[i*nCols + j] = j;
      }
    }
    rowOffsets[nRows] = nRows*nCols;
    return CSRMatrix(std::move(vals), std::move(rowOffsets), std::move(colInds), nCols);
  }

  template <class U> CSRMatrix<U>
  clone() const {
    CSRMatrix<U> cln(rowOffsets_, columns_);
    return cln;
  }

  template <class Lambda>
  CSRMatrix mutate(Lambda&& fxn) const {
    CSRMatrix m(rowOffsets_, columns_, nCols_);
    auto nelems = values_.size();
    for (size_t idx=0; idx < nelems; ++idx){
      m.values_[idx] = fxn(values_[idx]);
    }
    return m;
  }

  int colAt(int idx) const {
    return columns_[idx];
  }

  T& valueAt(int idx) {
    return values_[idx];
  }

  const T& valueAt(int idx) const {
    return values_[idx];
  }


  int rowOffset(int row) const {
    return rowOffsets_[row];
  }

  const int* colData() const {
    return columns_.data();
  }

  std::tuple<int,const T*, const int*> getRow(int row) const {
    int nCols = rowOffsets_[row+1] - rowOffsets_[row];
    int offset = rowOffsets_[row];
    return std::make_tuple(nCols, &values_[offset], &columns_[offset]);
  }

  std::tuple<int,T*,int*> getRow(int row) {
    int nCols = rowOffsets_[row+1] - rowOffsets_[row];
    int offset = rowOffsets_[row];
    return std::make_tuple(nCols, &values_[offset], &columns_[offset]);
  }

 private:
  CSRMatrix(const std::vector<int>& rows, const std::vector<int>& cols, int ncols) :
    values_(cols.size()), rowOffsets_(rows), columns_(cols), nCols_(ncols)
  {
  }

  std::vector<T> values_;
  std::vector<int> rowOffsets_;
  std::vector<int> columns_;
  int nCols_;

};

class HungarianAlgorithm
{
public:
  HungarianAlgorithm(const CSRMatrix<double>& m);

  double solve(const CSRMatrix<double>& m);

  const std::vector<int>& assignment() const {
    return coveredRows_;
  }

private:
  void findOptimal(CSRMatrix<double>& m);
  void findKnownOptimal(CSRMatrix<double>& m);
  double computeCost(const CSRMatrix<double>& m);
  void step2(CSRMatrix<double>& m);
  void step3(CSRMatrix<double>& m);
  void step4(CSRMatrix<double>& m);
  void visitRowStep3(CSRMatrix<double>& matrix, int row, std::deque<int>& rowsToVisit,
                     int& numCoveredRows, int& numCoveredCols);

 private:
  std::vector<bool> markedCols_;
  std::vector<bool> markedRows_;
  std::vector<int> coveredColumns_;
  //stores RELATIVE column numbers, ot absoluate column numbers
  std::vector<int> coveredRows_;
};


#endif
