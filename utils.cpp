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

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/hypergeometric.hpp>

#include "utils.h"
using namespace std;
using namespace boost::math;
using boost::math::chi_squared;
using boost::math::quantile;
using boost::math::complement;
using boost::math::cdf;

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
	cerr << "Execution time for operation : " << operation << " : "  << end.tv_usec - begin.tv_usec << " µs"  << endl;
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

double chisquare( vector<double> toTest , vector<double> all )
{
	boost::math::chi_squared chi(1);
	double a1 = toTest[0] ;
	double a2 = toTest[1];
	double b1 = all[0];
	double b2 = all[1];;


	double s = a1 + a2 + b1 + b2;
	double K = s * (a1 * b2 - a2 * b1) * (a1 * b2 - a2 * b1) / (a1 + a2) / (b1 + b2) / (a1 + b1) / (a2 + b2);
	double P = boost::math::cdf(chi, K);

	return P;
}


double fisher_test(vector<double> toTest, vector<double> control )
{
	double a = toTest[0];
	double b = toTest[1];
	double c = control[0];
	double d = control[1];

	double N = a + b + c + d;
	double r = a + c;
	double n = c + d;
	double max_for_k = min(r, n);
	double min_for_k = (double)max(0, int(r + n - N));
	hypergeometric_distribution<> hgd(r, n, N);
	double cutoff = pdf(hgd, c);
	double tmp_p = 0.0;
	for(int k = min_for_k;k < max_for_k + 1;k++)
	{
		double p = pdf(hgd, k);
		if(p <= cutoff)
		{
			tmp_p += p;
		}
	}
	return tmp_p;
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
	double sumone;
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

// Method for calculating fisher exact test 2-sided, return the pvalue.
double FET( int a , int b , int c , int d )
{
	int n = a + b + c + d;
	double logpCutOff = logHypergeometricProb( a , b , c , d );
	double pFraction = 0;
	double logpValue = 0;

	for( int x = 0 ; x <= n ; x ++ )
	{
		if( ( a + b - x >= 0 ) && ( a + c - x >= 0 ) && ( d - a + x >= 0 ) )
		{
			double l = logHypergeometricProb( x , a + b - x , a + c - x , d - a + x );
			if( l <= logpCutOff )
			{
				pFraction += exp( l - logpCutOff );
			}
		}
	}
	logpValue = logpCutOff + log( pFraction );

	return exp(logpValue);
}

// method for calculating the hypergeometrical log value  for the FET.
double logHypergeometricProb( int a , int b , int c  , int d )
{
	return logFactoriel( a + b ) + logFactoriel( c + d ) + logFactoriel( a + c ) + logFactoriel( b + d )- logFactoriel( a ) - logFactoriel( b ) - logFactoriel( c ) - logFactoriel( d ) - logFactoriel( a + b + c + d );
}

// Method for calculating a log factoriel
double logFactoriel( int inc )
{
	double ret;
	for( ret = 0 ; inc > 0 ; inc -- )
	{
		ret += log( (double)inc );
	}
	return ret;
}
