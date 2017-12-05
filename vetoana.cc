/*
* N. Warrack
* Program to analyse output from rivet analysis with letter codes for different veto reasons
* use, for example: 
* river -a ATLAS_2014_I1304289 fifo.hepmc >> vetoCodes.log
* where output looks like this (might need to use cout in your rivet analysis .cc file though!!):
* 1
* 1
* 1
* 3
* 1 
* 1
* 5
* 1
* etc.
* etc.
*/

#include <iostream>
#include <fstream>
#include "datafileprocess.h"

using namespace std;


int main(){

ifstream mydata;
  // Now open an existing file...
 mydata.open("ana.log");
 // check that operation worked...
 if ( mydata.fail() ) {
   cout << "Sorry, couldn’t open file ana.log -- It does not exisit in this location!" << endl;
   return 1; 
 }
 
 // Count lines in vetoCodes.log
 int lines = lineCounter(mydata);
 
 int vetoFrequency[] = {0,0,0,0,0,0,0,0,0} ;

 int x = 0 ;
 int line = 0 ;
 string infoline ;
 while ( mydata.good()){
     // Ignore first 20 lines to hopefully skip fastjets initialisation blurb (check vetoCodes.log by eye)
     if ( line < 30 ) {
       if( mydata.eof() ) break ;             
       getline(mydata, infoline);
       cout << "WARNING: ignored line -> " << infoline << endl;
       line++ ;
     } else {
       line++;
       mydata >> x ;
       //       cout << "x=" << x << endl;
          
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
       
     }
 }
 
   mydata.close(); // close when finished   
   cout << line-20 << " lines read from file" << endl;
   

   
   cout << "veto code 1: " << vetoFrequency[0] <<endl ;
   cout << "veto code 2: " << vetoFrequency[1] <<endl ;
   cout << "veto code 3: " << vetoFrequency[2] <<endl ;
   cout << "veto code 4: " << vetoFrequency[3] <<endl ;
   cout << "veto code 5: " << vetoFrequency[4] <<endl ;
   cout << "veto code 6: " << vetoFrequency[5] <<endl ;
   cout << "veto code 7: " << vetoFrequency[6] <<endl ;
   cout << "veto code 8: " << vetoFrequency[7] <<endl ;
   cout << "passed: " << vetoFrequency[8] <<endl ;
   
   return 0;
   
   }
