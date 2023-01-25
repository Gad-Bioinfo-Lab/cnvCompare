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
	bool useVCFFormat = true; 
	bool useBEDFormat = false;

	int currentThread = 0;
	int depthThreshold;
	int loggingLevel;
	int filterSize = 0;

	double ratioThreshold;
	double pvalueThreshold;


	string inputFile;
	string inputControlFile;
	string logFile;


	po::options_description desc("Allowed options");
	desc.add_options()
	("help,h", "produce help message")
	("input,i",  po::value<string>( &inputFile ), "List of input file(s) containing detected CNV from samples" )
	("control,c",  po::value<string>( &inputControlFile ), "List of input file(s) containing detected CNV from control" )
	("filter,f", po::value<int>( &filterSize ), "Minimum size for a CNV to be counted (0)" )
    ("whole,w" , "Whole mode. WARNING : Needs large amount of RAM" )
	("vcf" , "VCF mode : input files are in VCF format according to vcf specification v4.7 (default)")
	("bed", "BED mode : input files are in bed format + fields for cnv level and quality scores");


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

	cnvCompare * App;

	// check provided files
	if( inputControlFile.length( ) > 0 )
	{
		if( !IsFileReadable( inputControlFile ) )
		{
			cerr << "File provided as input : " << inputControlFile << " is not accessible : stopping" << endl;
			return -1;
		}
		App = new cnvCompare(inputFile, inputControlFile, 1, filterSize);
	}
	else
	{
		App = new cnvCompare(inputFile, 1, filterSize);
		cerr << "No file provided as input control file" << endl;
	}

	
  if (vm.count("vcf")) {
	useVCFFormat = true;
	useBEDFormat = false; 
  }
  if (vm.count("bed")) {
	useVCFFormat = false; 
	useBEDFormat = true;
  }

  App->setFormat(useVCFFormat , useBEDFormat);
	

  if( vm.count( "whole" ) ) {
    App->altLoop();
  } else {
    App->mainLoop();
  }



	delete App;
	cerr << "end of main" << endl;
  return 0;
}
