#include <iostream>
#include <vector>
#include <math.h>
using std::cout;
using std::cin;
using std::string;
using std::vector;

// lprogram to calculate the distance matrix from an input
// list of coordinates. first line of input file must have number
// of coordinates. must be 3 dimensional, each coordinate triple
// on one line.
// output will have number of lines equal to number of coordinates
// each line with number of coordinates values (a square array)

//////////////////////////////////////////////////////////////////////////////////////////
// functions 
//////////////////////////////////////////////////////////////////////////////////////////

// returns the magnitude of a vector
long double magnitude( const vector<long double> &vec ) {
    long double mag;
    
    mag = 0.0;
    for ( long double component : vec ) {
        mag += component*component;
    }
    return fsqrt( mag );
}
// pipe a vector to standard output
void coutVector( const vector<long double> &vec ) {
    for ( long double component : vec ) {    // use range ability for for-loop
    	cout << component << '\t';
    }
    cout << '\n';
}
// pipe an array to standard output
void coutArray( const vector<vector<long double>> &array ) {
    for ( vector<long double> vec : array ) {    // use range ability for for-loop
    	for ( long double component : vec ) {
	    	cout << component << '\t';
	}
	cout << '\n';
    }
}
// takes difference of two vectors
void diffVector( vector<long double> &vecResult, const vector<long double> &vec1, const vector<long double> &vec2 ) {
    for ( int i=0; i<vec1.size(); i++ ) {
        vecResult.at(i) = vec1.at(i) - vec2.at(i);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////
///// main progam
//////////////////////////////////////////////////////////////////////////////////////////
int main() {

     // declare variables with default/initial values
    int numCoordinates=0, dims=3;
    vector<long double> difference(dims,0.0);
  
    /////////////////////////////////////////////////////////////
    // parameter input from standard input
    /////////////////////////////////////////////////////////////
    
    // first read in the number of coordinates
    cin >> numCoordinates;
      
    // now declare data arrays with correct sizes
    vector<vector<long double>> x( numCoordinates, vector<long double>( dims, 0.0 ) );
    vector<vector<long double>> distanceMatrix( numCoordinates, vector<long double>( numCoordinates, 0.0 ) );
    
    for ( int i=0; i<numCoordinates; i++ ) {
    	for ( int j=0; j<dims; j++ ) {
    	    cin >> x.at(i).at(j);
        }
    }
     

//    cout << "# numCoordinates= " << numCoordinates << '\n';
//   cout << "# x = " ; coutArray(x); 
   
    /////////////////////////////////////////////////////////////
    // main algorithm
    ///////////////////////////////////////////////////////////// 
    
    for ( int i=0; i<numCoordinates-1; i++ ) {
    	for ( int j=i+1; j<numCoordinates; j++ ) {
//    	    cout << i << ' ' << j << '\n' << ' ';
 //   	    coutVector( x.at(i) );
  //  	     coutVector( x.at(j) );
    	    diffVector( difference, x.at(i), x.at(j) );
   // 	     coutVector( difference );
    	    distanceMatrix.at(i).at(j)=magnitude( difference );
    	    distanceMatrix.at(j).at(i)=magnitude( difference );
    	}
    }
     	    
    coutArray( distanceMatrix );

}
   
