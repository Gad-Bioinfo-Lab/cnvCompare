// C++ std libs
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <sys/time.h>
#include <sys/types.h>

// C++ Boost
#include <boost/program_options.hpp>

// utils
#include "utils.h"

// Classes
#include "cnvCompare.h"


using namespace std;
namespace po = boost::program_options;

int main(int argc, char* argv[])
{
  cerr << "Starting Main" << endl;

	bool countStop = false;

	int currentThread = 0;
	int depthThreshold;
	int loggingLevel;

	double ratioThreshold;
	double pvalueThreshold;


	string inputFile;

	string logFile;


	po::options_description desc("Allowed options");
	desc.add_options()
	("help,h", "produce help message")
	("input,i",  po::value<string>( &inputFile ), "Input TSV file(s) containing detected CNV from samples" )
    ("whole,w" , "Whole mode. WARNING : Needs large amount of RAM" );


	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);



	if( argc <= 1 )
	{
		cerr << "Error while checking program arguments" << endl;
		cerr << desc << "\n";
		return 1;
	}


	// check provided files
	if( inputFile.length( ) > 0 )
	{
		if( !IsFileReadable( inputFile ) )
		{
			cerr << "File provided as input : " << inputFile << " is not accessible : stopping" << endl;
			return -1;
		}
	}
	else
	{
		cerr << "No file provided as input file : stopping" << endl;
		return -1;
	}

	cnvCompare * App = new cnvCompare( inputFile , 1 );

  if( vm.count( "whole" ) ) {
    App->altLoop();
  } else {
    App->mainLoop();
  }



	delete App;
	cerr << "end of main" << endl;
  return 0;
}
