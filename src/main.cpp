#define BOOST_LOG_DYN_LINK 1

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
#include <execinfo.h>
#include <signal.h>

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

// namespaces
using namespace std;
using namespace boost;
namespace po = boost::program_options;
namespace logging = boost::log;
namespace keywords = boost::log::keywords;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;

// logging initiation function 
void init_logging() {
	logging::add_file_log(keywords::file_name ="cnvCompare_%N.log" , keywords::format = "[%TimeStamp%]: %Message%");
	logging::core::get()->set_filter(logging::trivial::severity >= logging::trivial::info);
}

// signal handler to leave properly if seg fault, or interruption. 
void handler(int sig) {
  void *array[10];
  size_t size;
  // get void*'s for all entries on the stack
  size = backtrace(array, 10);
  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

int main(int argc, char* argv[])
{
	// catching signals
	signal(SIGSEGV | SIGINT | SIGTERM | SIGABRT | SIGFPE, handler);

	// logging start
	init_logging();
  	BOOST_LOG_TRIVIAL(trace) << "Starting Main" << endl;

	// decla
	bool useVCFFormat = true; 
	bool useBEDFormat = false;
	int filterSize = 0;
	string inputFile;
	string inputControlFile;
	string logFile;
	string dictFile; 
	string suffix;

	// Option management
	po::options_description desc("Allowed options");
	desc.add_options()
	("help,h", "Produce help message")
	("version,v", "Print version and exit")
	("input,i",  po::value<string>( &inputFile ), "List of input file(s) containing detected CNV from samples")
	("control,c",  po::value<string>( &inputControlFile ), "List of input file(s) containing detected CNV from control")
	("filter,f", po::value<int>( &filterSize ), "Minimum size for a CNV to be counted (0)")
    ("whole,w", "Whole mode. WARNING : Needs large amount of RAM")
	("fast,q", "Fast mode. WARNING : experimental")
	("dict,d", po::value<string>( &dictFile ), "Dictionnary used to populate the chromosome list")
	("suffix,s", po::value<string>( &suffix ), "Suffix to use for the output files (default : count")
	("vcf", "VCF mode : input files are in VCF format according to vcf specification v4.7 (default)")
	("bed", "BED mode : input files are in bed format + fields for cnv level and quality scores");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	// deal with basic options 
	if( argc <= 1 )
	{
		cerr << "Error while checking program arguments" << endl;
		cerr << desc << "\n";
		return 1;
	}
	if (vm.count("help")) {
		cerr << desc << "\n";
		return 0;
	}
	if (vm.count("version")) {
		cerr << "cnvCompare : comparing and counting CNV's detected by sequencing experiments" << endl; 
		cerr << "Version 1.5.0" << endl; 
		cerr << "WARNING : do not use in diagnostics, results are not guaranteed" << endl;
		cerr << "Author : <yannis.duffourd@u-bourgogne.fr> - INSERM U1231 GAD - CHU Dijon" << endl; 
		cerr << "Copyright (c) 2022 GAD Lab under the aGPL v3 License" << endl;
		return 0;
	}

	// check provided files
	if( inputFile.length( ) > 0 ) {
		if( !IsFileReadable( inputFile ) ) {
			cerr << "File provided as input : " << inputFile << " is not accessible : stopping" << endl;
			return 1;
		}
	} else {
		cerr << "No file provided as input file : stopping" << endl;
		return 1;
	}

	// Declare the App
	cnvCompare * App;

	// check provided files and contruct object according to available options
	// consider using a template ? Is it possible in a constructor ?
	if( inputControlFile.length( ) > 0 ) {
		if( !IsFileReadable( inputControlFile ) ) {
			cerr << "File provided as control input : " << inputControlFile << " is not accessible : stopping" << endl;
			return 1;
		}
		App = new cnvCompare(inputFile, inputControlFile, 1, filterSize);
	} else {
		App = new cnvCompare(inputFile, 1, filterSize);
		BOOST_LOG_TRIVIAL(warning) << "No file provided as input control file" << endl;
	}

	// deal with the dict option
	if (vm.count("dict")) {
		if( dictFile.length( ) > 0 ) {
			if( !IsFileReadable( dictFile ) ) {
				cerr << "File provided as dictionnary : " << dictFile << " is not accessible : using default parameters" << endl;
				App->setHasDict(false);
			} else { 
				App->setDictFile(dictFile); 
				App->setHasDict(true);
			}
		}
	} else {
		App->setHasDict(false);
	}

	// configure file format
	if (vm.count("vcf")) {
		useVCFFormat = true;
		useBEDFormat = false; 
	}
	if (vm.count("bed")) {
		useVCFFormat = false; 
		useBEDFormat = true;
	}
	App->setFormat(useVCFFormat , useBEDFormat);
	if (vm.count("suffix")) {
		App->setSuffix(suffix); 
	}

	// configure running mode and launch the appropriate loop
	if( vm.count( "whole" ) ) {
		App->altLoop();
	} else {
		if (vm.count("fast")) {
			App->fastLoop();
		} else {
			App->mainLoop();
		}
	}

	// GC
	delete App;
	BOOST_LOG_TRIVIAL(trace) << "end of main" << endl;
	return 0;
}
