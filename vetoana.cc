/*
* N. Warrack
* Program to analyse output from rivet analysis with letter codes for different veto reasons
* use, for example: 
* river -a ATLAS_2014_I1304289 fifo.hepmc > vetoCodes.log
* where digits followed by a space are cout'ed from a rivet analysis and form a line for each event
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>

using namespace std;


int lineCounter(ifstream& dat) {
  int ctr = 0 ;
  string line ;
  while (getline(dat, line)) ctr++ ;
  // clear tags and return ifstream to beginning of input file
  dat.clear();
  dat.seekg(0, ios::beg);
  return ctr ;
}


template<typename T>
void printout(T const &rhs) {
  cout << rhs << endl ;
}



int main ( int argc, char *argv[] ){

  // check for appropriate number of command line arguments  
  if (argc != 3) {
    cout << "ERROR: requires an input (.log) file and an integer (number of lines to skip) as arguments... " << endl ;
    cout << "ERROR: (see input file - count number of lines to skip by eye)" << endl ; 
    cout << "ERROR: ...aborting."<< endl;
    return 1;
  }


  // feedback user info:
  cout << "INFO: input file: " ; 
  printout(argv[1]) ;
  cout << "INFO: integer number of lines to skip: " ; 
  printout(argv[2]) ;


 // Open an existing file given as argument 
  string filename = argv[1] ; 
  ifstream mydata ;
  mydata.open(filename) ;
  if ( mydata.fail() ) {
    cout << "ERROR: couldn’t open file: " << filename << endl ;
    return 1; 
  }

  
  // Count lines in input file
  int const static lines = lineCounter(mydata) ;
  int vetoFrequency[] = {0,0,0,0,0,0,0,0,0} ;
  const static int m = atoi(argv[2]) ;
  int x = 0 ;
  int line = 0 ;
  string linestring ;
  

  // abort of lines to skip > lines in file
  if (m > lines) {
   cout << "ERROR: line to skip is larger than lines in file " << endl ; 
   return 1;
  }


  // Begin file analysis
  while ( mydata.good()){
    
    if ( line < m ) {
      if( mydata.eof() ) break ;             
      getline(mydata, linestring);
      // cout << "WARNING: ignoring first " << m << " lines..." << line + 1 << " -> " << infoline << endl; // uncomment this line to print skipped lines to screen
      line++ ;
      
    } else { // if ( line >= m )
      
      // output info to user
      cout << "INFO: Ignored first " << m << " lines..." << endl ;
      cout << "INFO: reading remaining " << lines - m << " lines from file" << endl ;
      
      while(getline(mydata, linestring))
	{
	  istringstream iss(linestring) ;
	  while(iss >> x){
	    //___________________________________________________________________________
	    if ( x == 1 ) {vetoFrequency[0] =(vetoFrequency[0]+1);}else{
	      if ( x == 2 ) {vetoFrequency[1] =(vetoFrequency[1]+1);}else{
		if ( x == 3 ) {vetoFrequency[2] =(vetoFrequency[2]+1);}else{
		  if ( x == 4 ) {vetoFrequency[3] =(vetoFrequency[3]+1);}else{
		    if ( x == 5 ) {vetoFrequency[4] =(vetoFrequency[4]+1);}else{
		      if ( x == 6 ) {vetoFrequency[5] =(vetoFrequency[5]+1);}else{
			if ( x == 7 ) {vetoFrequency[6] =(vetoFrequency[6]+1);}else{
			  if ( x == 8 ) {vetoFrequency[7] =(vetoFrequency[7]+1);}else{
			    if ( x == 9 ) {vetoFrequency[8] =(vetoFrequency[8]+1);}
			   }
			 }
		       }
		     }
		   }
		 }
	       }
	     }
	     //___________________________________________________________________________
	   } // end of while statement
	   
	   line++;
	 }
     }
 }
 
 mydata.close(); // close when finished   
 
 cout << "veto code 1: " << vetoFrequency[0] <<endl ;
 cout << "veto code 2: " << vetoFrequency[1] <<endl ;
 cout << "veto code 3: " << vetoFrequency[2] <<endl ;
 cout << "veto code 4: " << vetoFrequency[3] <<endl ;
 cout << "veto code 5: " << vetoFrequency[4] <<endl ;
 cout << "veto code 6: " << vetoFrequency[5] <<endl ;
 cout << "veto code 7: " << vetoFrequency[6] <<endl ;
 cout << "veto code 8: " << vetoFrequency[7] <<endl ;
 cout << "PASSED (=9): " << vetoFrequency[8] <<endl <<endl;
 
 return 0;
 
}
