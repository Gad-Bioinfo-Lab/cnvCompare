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
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>

// utils
#include "utils.h"

// Classes
#include "cnvCompare.h"


using namespace std;
using namespace boost;
namespace po = boost::program_options;
namespace logging = boost::log;
namespace keywords = boost::log::keywords;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;


void init_logging() {
	logging::add_file_log(keywords::file_name ="cnvCompare_%N.log" , keywords::format = "[%TimeStamp%]: %Message%");
	logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::debug);
	// cerr << "End of init logging" << logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::trace) << endl; 
}


int main(int argc, char* argv[])
{
	init_logging();
  	BOOST_LOG_TRIVIAL(trace) << "Starting Main" << endl;

	bool useVCFFormat = true; 
	bool useBEDFormat = false;

	// int currentThread = 0;
	int filterSize = 0;

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
			return 1;
		}
	}
	else
	{
		cerr << "No file provided as input file : stopping" << endl;
		return 1;
	}

	cnvCompare * App;

	// check provided files
	if( inputControlFile.length( ) > 0 )
	{
		if( !IsFileReadable( inputControlFile ) )
		{
			cerr << "File provided as input : " << inputControlFile << " is not accessible : stopping" << endl;
			return 1;
		}
		App = new cnvCompare(inputFile, inputControlFile, 1, filterSize);
	}
	else
	{
		App = new cnvCompare(inputFile, 1, filterSize);
		BOOST_LOG_TRIVIAL(warning) << "No file provided as input control file" << endl;
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
	BOOST_LOG_TRIVIAL(trace) << "end of main" << endl;
  return 0;
}
