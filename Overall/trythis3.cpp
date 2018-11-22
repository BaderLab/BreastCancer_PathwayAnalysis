/*
  This code is intended to be invoked by R code.  When attempting to
  run this on a GPC compute node, be sure to load the "extras" module,
  or else the 'make' command will not be available, and the runtime
  compilation of this code will fail.
*/

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix trythis3(IntegerVector data0, IntegerVector status0) {

// get the size of the data.
  int n = data0.size();

// the output matrix.
  IntegerMatrix M(2,2);

  int masterTable1 = 0;
  int caseSum = 0;
  int controlSum = 0;
  int caseSumMin = 0;
  int controlSumMin = 0;
  int tempSum = 0;
  
// do the sums.
  for(int i = 0; i < n; i++) {
    
    if (not ((IntegerVector::is_na(data0[i])) ||
	     (IntegerVector::is_na(status0[i])))) {
      
      
      caseSum += status0[i] * data0[i];
      controlSum += ((1-status0[i]) * data0[i]);
      tempSum += (1-status0[i]);
	  masterTable1 += status0[i];
	  
    }

  }

// create the values.
  caseSumMin = masterTable1 * 2 - caseSum;
  controlSumMin = tempSum * 2 - controlSum;

// populate the output.
  M(0,0) = controlSum;
  M(0,1) = caseSum;
  M(1,0) = controlSumMin;
  M(1,1) = caseSumMin;

  return M;

}
