// utility functions for bioinformatics
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <sys/time.h>
#include <sys/types.h>

#include "utils.h"
using namespace std;


// test if a file is readable
// return value : true if readable ; false if not
bool IsFileReadable( string file )
{
	ifstream fichier( file.c_str() );
	return !fichier.fail();
}

// display a time lenth in µs
// return value : void
void ExecMeasure( struct timeval begin , struct timeval end , string operation )
{
	
	string unit = "µs";
	if ((end.tv_usec - begin.tv_usec) < 0) {
		end.tv_usec /= 1000; 
		begin.tv_usec /= 1000;
		unit = "ms";
	}
	// if ((end.tv_usec - begin.tv_usec) > sizeof(suseconds_t)) {
	// 	end.tv_usec /= 1000; 
	// 	begin.tv_usec /= 1000;
	// 	unit = "s";
	// }
	cerr << "Execution time for operation : " << operation << " : "  << end.tv_usec - begin.tv_usec << " " << unit << endl;
}


int string_to_int( string incomingStr )
{
	istringstream isstmp( incomingStr );
	int i;
	isstmp >> i;
	return i;
}

string double_to_string( double incoming )
{
	string result;
	ostringstream oss;
	oss << incoming;
	result = oss.str();
	return result;
}


string int_to_string( int incoming )
{
	string result;
	ostringstream oss;
	oss << incoming;
	result = oss.str();
	return result;
}

string pyReplace( string incoming , string pattern , string replacement )
{
	while (incoming.rfind( pattern ) != string::npos )
	{
		int n = incoming.rfind( pattern );
		int l = pattern.length();

		incoming.replace( n , l , replacement );
	}
	return incoming;
}

string char_to_string(char incoming)
{
	string s;
	stringstream ss;
	ss << incoming;
	ss >> s;
	return s;
}

vector<string> parseOnSep( string inc , string sep )
{
	// cerr << "Entering ParseOnSep function" << endl;
	// cerr << "\tIncoming string : " << inc << " ; separator : " << sep << endl;

	vector<string> ret;
	istringstream issInc (inc);
	string mot;
	while ( getline( issInc, mot, string_to_char( sep ) ) )
	{
		ret.push_back( mot );
	}
	return ret;

}


char string_to_char( string inc )
{
	char cstr[inc.size() + 1];
	inc.copy(cstr, inc.size() + 1);
	cstr[inc.size()] = '\0';
	return *cstr;
}




string strip( string inc )
{
	cerr << "Passing into strip << " << inc ;
	string::size_type pos = 0;
	while ( ( pos = inc.find ("\n",pos) ) != string::npos )
	{
		cerr << " ; pos = " << pos ;
		inc.erase ( pos, 2 );
	}
	cerr << " to " << inc << endl;
	return inc;
}




char checkBase(char incoming )
{
	if( incoming == 'c' )
	{
		return 'C';
	}
	if( incoming == 't' )
	{
		return 'T';
	}
	if( incoming == 'a' )
	{
		return 'A';
	}
	if( incoming == 'g' )
	{
		return 'G';
	}
	if( incoming == 'n' )
	{
		return 'N';
	}
	if( incoming == 'C' )
	{
		return 'C';
	}
	if( incoming == 'T' )
	{
		return 'T';
	}
	if( incoming == 'A' )
	{
		return 'A';
	}
	if( incoming == 'G' )
	{
		return 'G';
	}
	if( incoming == 'N' )
	{
		return 'N';
	}
	return 'N';
}



//Method for calculating a sd from a vector of double
double sd_calculator( vector<double> incVector )
{
	// Déclarations
	double sd;
	double temp_value;
	double sumone = 0;
	double sumtwo = 0;
	double moyenne;
	int number = 0;
	double variance;

	// calcul des moyennes et moyennes carrées
	vector<double>::iterator myIter;
	for( myIter = incVector.begin() ; myIter != incVector.end() ; myIter ++ )
	{
		temp_value = * myIter;

		sumone += temp_value;
		sumtwo += ( temp_value * temp_value );
		number ++;
	}

	// calcul de la moyenne
	moyenne = sumone / number;
	// Calcul de la variance
	variance = ( sumtwo / number ) - ( moyenne * moyenne );
	// Calcul ecart type
	sd = sqrt(variance);

	return sd;
}

// Method for calculating a mean from a vector of double
double moyenne_calculator( vector<double> incVector )
{
	// Déclarations
	double temp_value;
	double sumone = 0.0;
	double moyenne;
	int number = 0;

	// calcul des moyennes et moyennes carrées
	vector<double>::iterator myIter;
	for( myIter = incVector.begin() ; myIter != incVector.end() ; myIter ++ )
	{
		temp_value = * myIter;
		sumone += temp_value;
		number ++;
	}

	// calcul de la moyenne
	if( number != 0 )
	{
		moyenne = sumone / number;
	}
	else
	{
		return 0;
	}
	return moyenne;
}

