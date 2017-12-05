#ifndef DATAFILEPROCESS_H
#define DATAFILEPROCESS_H

using namespace std;

// process to count dimensionality of input data file
int lineCounter(ifstream& dat){

  int ctr = 0 ;
  string line ;
  while (getline(dat, line)) ctr++ ;

  /// return ifstream to beginning of input file
  dat.clear();
  dat.seekg(0, ios::beg);

  return ctr ;
}

#endif
