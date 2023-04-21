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

double sd_calculator( std::vector<double> );
double moyenne_calculator( std::vector<double> );
double moyenne_calculator_array( double , int );
#endif
