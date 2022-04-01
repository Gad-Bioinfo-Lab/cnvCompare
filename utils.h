#ifndef _UTILS_H
#define _UTILS_H

#include <vector>


bool IsFileReadable( std::string );
void ExecMeasure( struct timeval, struct timeval, std::string );
int string_to_int( std::string );
std::string double_to_string( double );
std::string int_to_string( int );
char checkBase(char);
std::string pyReplace( std::string , std::string , std::string );
std::string char_to_string(char);
std::vector<std::string> parseOnSep( std::string , std::string );
char string_to_char( std::string );
std::string strip( std::string );
double chisquare( std::vector<double>, std::vector<double> );
double fisher_test(std::vector<double> , std::vector<double> );
double sd_calculator( std::vector<double> );
double moyenne_calculator( std::vector<double> );
double FET( int , int , int , int );
double logHypergeometricProb( int , int , int , int );
double logFactoriel( int );

#endif
