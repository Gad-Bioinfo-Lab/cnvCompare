#define BOOST_LOG_DYN_LINK 1

// C++ std libs
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/time.h>
#include <sys/types.h>
#include <vector>
#include <list>
#include <iomanip>
#include <unordered_map>
#include <filesystem>
//#include <ranges>

// Boost 
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/expressions.hpp>
#include <boost/filesystem/path.hpp>

// utils
#include "utils.h"

// Classes
#include "cnvCompare.h"

// namespaces
using namespace std;
using namespace boost;
namespace logging = boost::log;
namespace fs = boost::filesystem;

/**
 * @brief default constructor (useless)
 * @param none
 * @return none
 **/
cnvCompare::cnvCompare() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::cnvCompare default constructor" << endl;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::cnvCompare default constructor" << endl;
}


/**
 * @brief Constructor of the class without controls
 * @param iF file containing the list of input files
 * @param iF number of threads 
 * @param s size of the filter to apply
 * @return none
 **/
cnvCompare::cnvCompare(string iF, int nT, int s) {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::cnvCompare (string, int, int) Ctor" << endl;
  BOOST_LOG_TRIVIAL(trace) <<  "cnvCompare constructor called";
  BOOST_LOG_TRIVIAL(trace) <<  "sizeof(char) = " << sizeof(char);
  BOOST_LOG_TRIVIAL(trace) <<  "sizeof(string) = " << sizeof(string);
  BOOST_LOG_TRIVIAL(trace) <<  "sizeof(int32_t) = " << sizeof(int32_t);
  BOOST_LOG_TRIVIAL(trace) <<  "sizeof(int) = " << sizeof(int);
  
  // init of the needed class members
  this->inputFile = iF;
  this->nbThread = nT;
  this->useControls = false;
  this->filterSize = s;
  this->suffix = "count";
  // Adding sample files 
  int n = this->fillMap(iF, "sample");
  BOOST_LOG_TRIVIAL(info) << n << " files added to sample list" << endl;
  // populating chromosome either by parsing the dict of with a default human map
  if (this->hasDict) {
    this->parseDictFile(this->getDictFile());
    
  } else {
    int ret = this->populateChr();
    if (ret != 0) {
      BOOST_LOG_TRIVIAL(error) << "Chromosome population couldn't be filled in" << endl;
    }
  }
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::cnvCompare (string, int, int) Ctor" << endl; 
}

/**
 * @brief Constructor of the class with controls
 * @param iF file containing the list of input files
 * @param cF file containing the list of control files
 * @param iF number of threads 
 * @param s size of the filter to apply
 * @return none
 **/
cnvCompare::cnvCompare(string iF, string cF, int nT, int s) {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::cnvCompare (string, string, int, int) Ctor" << endl;
  // init of the needed class members
  this->inputFile = iF;
  this->controlFile = cF;
  this->nbThread = nT;
  this->useControls = true;
  this->filterSize = s;
  this->suffix = "count";
  // Adding sample files & control files
  int n = this->fillMap(iF, "sample");
  BOOST_LOG_TRIVIAL(info) << n << " files added to sample list" << endl;
  int m = this->fillMap(cF, "control");
  BOOST_LOG_TRIVIAL(info) << m << " files added to control list" << endl;
  // populating chromosome either by parsing the dict of with a default human map
  if (this->hasDict) {
    this->parseDictFile(this->getDictFile());
  } else {
    int ret = this->populateChr();
    if (ret != 0) {
      BOOST_LOG_TRIVIAL(error) << "Chromosome population couldn't be filled in" << endl;
    }
  }
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::cnvCompare (string, string, int, int) Ctor" << endl;
}


/**
 * @brief Main loop used to get data, and compute counts chr by chr
 * @param none
 * @return none
 **/
void cnvCompare::mainLoop() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::mainLoop " << endl;
  for (auto &a : this->chromosomeMap) {
    // getting data
    this->getDatabyChr(a);
    // computing counts 
    this->computeCountsbyChr(a);
    // GC 
    this->cleanData();
  }
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::mainLoop " << endl;
}

/**
 * @brief Alt loop used to get data, and compute counts on the whole genome ; needs huge amount of RAM
 *        Warning : this method is deprecated. Do not use. 
 * @param none
 * @return none
 **/
void cnvCompare::altLoop() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::altLoop " << endl;
  // getting data
  this->getDataWhole();
  // computing counts 
  this->computeCountsWhole();
  // GC 
  this->cleanData();
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::altLoop " << endl;
}

/**
 * @brief fast loop used to get data, and compute counts on intervals ; experimental feature
 * @param none
 * @return none
 **/
void cnvCompare::fastLoop() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::fastLoop " << endl;
  // getting data
  this->getDataFast();
  // computing counts 
  this->computeCountsFast();
  // GC 
  this->cleanData();
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::fastLoop " << endl;
}


/**
 * @brief Method used to emptying data for the chromosomes
 * @param none
 * @return none
 **/
void cnvCompare::cleanData() { 
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::cleanData " << endl;
  this->dataByChr.clear();
  this->data.clear();
  this->breakpoints.clear(); 
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::cleanData " << endl; 
}

/**
 * @brief Method used to collect data from input files
 * @param incChr the chr to filter
 * @return none
 **/
void cnvCompare::getDatabyChr(string incChr) {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getDatabyChr " << endl;
  BOOST_LOG_TRIVIAL(info) << "Gathering data for chr " << incChr << endl;
  string ligne;
  string ligneCNV;
  string mot;
  string header = "#";
  string currentChr;
  string chromosome;
  string s_type;
  string s_start;
  string s_end;
  string s_value;

  // prefill the container
  unordered_map<long, short> tempMap2;
  this->dataByChr[0] = tempMap2;
  this->dataByChr[1] = tempMap2;
  this->dataByChr[2] = tempMap2;
  this->dataByChr[3] = tempMap2;
  this->dataByChr[4] = tempMap2;
  this->dataByChr[5] = tempMap2;

  this->nbFile = 0;

  // tsv parsing
  map<string, string>::iterator myIterA;
  for (myIterA = this->fileMap.begin(); myIterA != this->fileMap.end(); myIterA++) {
    ligne = myIterA->first;
    ifstream cnvStream(ligne.c_str());
    BOOST_LOG_TRIVIAL(info) << "\tReading file " << ligne << "\t" << endl;
    this->nbFile++;
    long nbLigneFile = 0;
    short nbOfConcernedIndividual = 0;
    while (getline(cnvStream, ligneCNV)) {

      // need to deal with header "#"
      if (ligneCNV.find(header) == 0) {
        continue;
      }

      vector<string> res;
      if (this->getFormat() == "BED")
      {
        res = this->parseBEDLine(ligneCNV);
      }
      else
      {
        res = this->parseVCFLine(ligneCNV);
        nbOfConcernedIndividual = string_to_int(res[5]);
      }
      
      // type conversion
      chromosome = res[0];

      // Pass the line if not the asked chromosome
      if (chromosome != incChr)
      {
        continue;
      }

      s_type = res[3];
      // Pass if not del / dup
      if ((s_type != "DEL") && (s_type != "DUP"))
      {
        continue;
      }
      long start = string_to_int(res[1]);
      long end = string_to_int(res[2]);
      unsigned int value = string_to_int(res[4]);

      // size filter
      if ((end - start) < this->getFilterSize())
      {
        continue;
      }

      // counting events after all filters
      nbLigneFile++;

      // thresholding the dups
      if (value > 5)
      {
        value = 5;
      }

      // inserting data 
      for (long i = start; i <= end; i++)
      {
        int value_to_insert = 0; 
        if ((this->dataByChr[value]).find(i) != (this->dataByChr[value]).end())
        {
          value_to_insert += this->dataByChr[value][i] + 1; 
        }
        this->dataByChr[value].insert_or_assign(i, value_to_insert);
      }
    }
    BOOST_LOG_TRIVIAL(info) << nbLigneFile << " events detected " << endl;
  }
  BOOST_LOG_TRIVIAL(info) << "Ended with " << this->getNbFile() << " files" << endl;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getDatabyChr " << endl;
}


/**
 * @brief Method used to compute counts from in memory data chr by chr
 * @param incChr the chr to filter
 * @return none
 **/
void cnvCompare::computeCountsbyChr(string incChr) {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::computeCountsbyChr " << endl;
  std::cout.precision(3);
  BOOST_LOG_TRIVIAL(info) << "Computing counts for chr " << incChr << endl;
  // struct timeval tbegin, tend;
  string ligne;
  string ligneCNV;
  string mot;
  string header = "#";
  int i = 0;
  string currentChr;
  string chromosome;
  string s_type;
  string s_start;
  string s_end;
  string s_value;
  string outFileName;

  // tsv parsing
  map<string, string>::iterator myIterA;
  for (myIterA = this->fileMap.begin(); myIterA != this->fileMap.end(); myIterA++) {
    // struct timeval tbegin, tend;
    ligne = myIterA->first;
    string s = myIterA->second;
    if (s == "control") {
      continue;
    }
    ifstream cnvStream(ligne.c_str());
    BOOST_LOG_TRIVIAL(info) << "\tReading file " << ligne << endl;

    // managing output file names
    fs::path pathObj(ligne);
    if (pathObj.has_extension()) {
      string extension = pathObj.extension().string(); 
      BOOST_LOG_TRIVIAL(debug) << "Extension detected : " << extension << endl;
      outFileName = pyReplace(ligne, extension, "." + this->getSuffix() + extension);
    } else {
      outFileName = ligne + "." + this->getSuffix();
    }
    if (outFileName == ligne) {
      BOOST_LOG_TRIVIAL(error) << "The output filename " << outFileName << " is the same as the input " << ligne << " : it will replace the original file : Stopping execution" << endl;
      exit(1);
    }
    BOOST_LOG_TRIVIAL(info) << "\t\tWriting in file " << outFileName << " in format : " << this->getFormat() << endl;
    // opening file in the apporpriate mode
    ofstream outStream;
    if (incChr == "chr1") {
      outStream.open(outFileName.c_str(), ios::out);
    } else {
      outStream.open(outFileName.c_str(), ios::out | ios::app);
    }

    // parsing input file
    while (getline(cnvStream, ligneCNV)) {
      // need to deal with header "#"
      if (ligneCNV.find(header) == 0) {
        if ((this->getFormat() == "VCF") && (incChr == "chr1")) {
          outStream << ligneCNV << endl;
        }
        continue;
      }
      vector<string> res;
      if (this->getFormat() == "BED") {
        res = this->parseBEDLine(ligneCNV);
      } else {
        res = this->parseVCFLine(ligneCNV);
      }

      // type conversion
      chromosome = res[0];
      s_type = res[3];
      long start = string_to_int(res[1]);
      long end = string_to_int(res[2]);
      unsigned int value = string_to_int(res[4]);
      if (value > 5) {
        value = 5;
      }
      // skip if not the good chr
      if (chromosome != incChr) {
        continue;
      }

      // pass if not del or dup
      BOOST_LOG_TRIVIAL(trace) << "s_type = " << s_type << endl;
      if ((s_type != "DUP") && (s_type != "DEL")) {
        outStream << ligneCNV << endl;
        continue;
      }

      // counts 
      double total = 0;
      for (long i = start; i <= end; i++) {
        total += (double)(this->dataByChr[value][i]);
      }
      double mean = total / (double)((end-start)+1);

      // need to adapt the output according to the choosen format
      if (this->getFormat() == "BED") {
        outStream << chromosome << "\t" << start << "\t" << end << "\t";
        if (value > 2) {
          outStream << "DUP\t";
        } else {
          outStream << "DEL\t";
        }
        outStream << value << "\t" << mean << "/" << this->getNbFile() << endl;
      } else {
        // output VCF
        res = this->parseVCFLine(ligneCNV);
        istringstream issLigne(ligneCNV);
        istringstream issInfo;
        string mot;
        string info;
        string infomot;
        string svtype;
        string ciend;
        string value;
        i = 0;

        while (getline(issLigne, mot, '\t')) {
          switch (i) {
          case 0:
            outStream << mot;
            break;
          case 7:
            outStream << "\t";
            info = mot;
            issInfo.str(info);
            while (getline(issInfo, infomot, ';')) {
              if (infomot.find("SVTYPE=") == 0) {
                svtype = parseOnSep(infomot, "=")[1];
                continue;
              }
              if (infomot.find("END=") == 0) {
                ciend = parseOnSep(infomot, "=")[1];
                continue;
              }
              if (infomot.find("VALUE=") == 0) {
                value = parseOnSep(infomot, "=")[1];
                continue;
              }
              outStream << infomot << ";";
            }
            outStream << "END=" << ciend << ";VALUE=" << value << ";SVTYPE=";
            if (string_to_int(value) > 2) {
              outStream << "DUP;";
            } else {
              outStream << "DEL;";
            }
            outStream << "COUNT=" << mean << "/" << this->getNbFile();
            break;
          default:
            outStream << "\t" << mot;
            break;
          }
          i++;
        }
        outStream << endl;
      }
      // gettimeofday(&tend, NULL);
      // ExecMeasure(tbegin, tend, "VCF writing");
    }
    outStream.close();
  }
  this->dataByChr.clear();
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::computeCountsbyChr " << endl;
}


/**
 * @brief Method used to parse a line from a BED file
 * @param incLine the line to parse
 * @return vector<string> containing values : chr, start, end, type, value
 * @todo need to check that all the values are present
 **/
// output a vector containing : chr start end type value from a BED line
vector<string> cnvCompare::parseBEDLine(string incLine) {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::parseBEDLine " << endl;
  vector<string> output;
  string mot;
  short int i = 0;

  // get information from the line
  istringstream issLigne(incLine);
  while (getline(issLigne, mot, '\t')) {
    switch (i) {
    case 0:
      output.push_back(mot);
      break;
    case 1:
      output.push_back(mot);
      break;
    case 2:
      output.push_back(mot);
      break;
    case 3:
      output.push_back(mot);
      break;
    case 4:
      output.push_back(mot);
      break;
    default:
      break;
    }
    i++;
  }
  return output;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::parseBEDLine " << endl;
}


/**
 * @brief Method used to parse a line from a VCF file
 * @param incLine the line to parse
 * @return vector<string> containing values : chr, start, end, type, value
 * @todo need to check that all the values are present
 **/
vector<string> cnvCompare::parseVCFLine(string incLine) {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::parseVCFLine " << endl;
  vector<string> output;
  map<string, string> temp;
  string mot;
  string infomot;
  string info;
  short GTindex = -1;
  short CNindex = -1;
  short int i = 0;
  double CNValue_d = 0.0; 
  int CNValue_i = 0;
  bool valueFound = false; 
  bool passGT = false;
  int nbOfConcernedIndiv = 0;
  vector<short> counts(6, 0); 
  istringstream issInfo;
  temp["SVTYPE"] = "NONE";
  vector<string> GTInfo;

  // get information from the line
  istringstream issLigne(incLine);
  i = 0;
  while (getline(issLigne, mot, '\t')) {
    if (i == 0) {
      output.push_back(mot);
      ++ i;
      continue;
    }
    if (i == 1) {
      output.push_back(mot);
      ++i;
      continue;
    }
    if (i == 7) {
      info = mot;
      issInfo.str(info);
      while (getline(issInfo, infomot, ';')) {
        if (infomot.find("SVTYPE=") == 0) {
          temp["SVTYPE"] = parseOnSep(infomot, "=")[1];
        }
        if (infomot.find("END=") == 0) {
          temp["END"] = parseOnSep(infomot, "=")[1];
        }
        if (infomot.find("VALUE=") == 0) {
          temp["VALUE"] = parseOnSep(infomot, "=")[1];
          valueFound = true;
        }
      }
      ++i; 
      continue;
    }
    if (i == 8) {
      // getting index for the GT field
      GTInfo = parseOnSep(mot, ":");
      for (long unsigned int n = 0 ; n <= GTInfo.size() ; n ++) {
        if (GTInfo[n] == "GT") {
          GTindex = n;
          BOOST_LOG_TRIVIAL(debug) << "GT index was found : " << GTindex << endl; 
        }
        if (! valueFound) {
          if (GTInfo[n] == "CN") {
            CNindex = n;
            BOOST_LOG_TRIVIAL(debug) << "CN index was found : " << GTindex << endl; 
          }
        }
      }
      ++i;
      continue;
    }
    
      // getting the non wild indiv 
    if (i >= 9) {
      // first checking if the GT field has been found
      if (passGT) {
        i++; 
        continue;
      }
      if ((temp["SVTYPE"] != "DEL") and (temp["SVTYPE"] != "DUP")) {
        break; 
      }

      if ((! valueFound) && (CNindex == -1)) {
        passGT = true;
        BOOST_LOG_TRIVIAL(info) << "No copy number value found on the VCF line " << incLine << " : passing it. Please check the VCF specifications" << endl;
        break;
      }
      if (GTindex == -1) {
        BOOST_LOG_TRIVIAL(info) << "No GT Found on the VCF line : counting only 1" << endl;
        nbOfConcernedIndiv = 1;
        passGT = true;
        break;
      } else {
        string GT = parseOnSep(mot, ":")[GTindex];
        BOOST_LOG_TRIVIAL(debug) << "\tGT Found : " << GT << endl;
        if ((GT != "0/0") and (GT != "0|0")) {
          // need to determine the copy level 
          if (! valueFound) {
            CNValue_d = stod(parseOnSep((parseOnSep(mot, ";")[CNindex]), "=")[1]);
            CNValue_i = floor(CNValue_d + 0.5);
            // not interested in cnv at n=2
            if (CNValue_i == 2) {
              ++i; 
              continue; 
            }
            if (CNValue_i > 5) {
              CNValue_i = 5;
            }
            counts[CNValue_i] += 1;
          }
          nbOfConcernedIndiv += 1;
          BOOST_LOG_TRIVIAL(debug) << "\t\tadding 1 concerned individual : " << GT << endl;
        }
      } 
      i++; 
      continue;
    }
    i++;
  }

  // add data to the vector from the temp map
  output.push_back(temp["END"]);
  output.push_back(temp["SVTYPE"]);
  if (! valueFound) {
    temp["VALUE"] = "-1";
  }
  output.push_back(temp["VALUE"]);
  output.push_back(int_to_string(nbOfConcernedIndiv));
  BOOST_LOG_TRIVIAL(debug) << "\tNumber of concerned individual is " << nbOfConcernedIndiv << endl;
  // transforming the counts vector into string 
  output.push_back(int_to_string(counts[0]) + "," + int_to_string(counts[1]) + "," + int_to_string(counts[2]) + "," + int_to_string(counts[3]) + "," + int_to_string(counts[4]) + "," + int_to_string(counts[5]));


  BOOST_LOG_TRIVIAL(debug) << "\tWill return output : " << endl;
  vector <string>::iterator myIter; 
  for (myIter = output.begin() ; myIter != output.end() ; myIter++ ) {
    BOOST_LOG_TRIVIAL(debug) << "\t\t" << *myIter << endl; 
  }
  BOOST_LOG_TRIVIAL(debug) << "\n";

  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::parseVCFLine " << endl;
  return output;
}

/**
 * @brief Method used to collect data from input files : Whole mode
 * @param none
 * @return none
 **/
void cnvCompare::getDataWhole() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getDataWhole " << endl;
  BOOST_LOG_TRIVIAL(info) << "Gathering data" << endl;
  // struct timeval tbegin, tend;
  string ligne;
  string ligneCNV;
  string mot;
  string header = "#";
  string currentChr;
  string chromosome;
  string s_type;
  string s_start;
  string s_end;
  string s_value;

  // tsv parsing
  map<string, string>::iterator myIterA;
  for (myIterA = this->fileMap.begin(); myIterA != this->fileMap.end(); myIterA++)
  {
    ligne = myIterA->first;
    ifstream cnvStream(ligne.c_str());
    BOOST_LOG_TRIVIAL(info) << "\tReading file " << ligne << "\t";
    this->nbFile++;
    long nbLigneFile = 0;
    while (getline(cnvStream, ligneCNV))
    {
      // need to deal with header "#"
      if (ligneCNV.find(header) == 0)
      {
        this->watchHeader(ligneCNV);
      }
      nbLigneFile++;
      vector<string> res;
      if (this->getFormat() == "BED")
      {
        res = this->parseBEDLine(ligneCNV);
      }
      else
      {
        res = this->parseVCFLine(ligneCNV);
      }

      // type conversion
      chromosome = res[0];
      s_type = res[3];
      long start = string_to_int(res[1]);
      long end = string_to_int(res[2]);
      unsigned int value = string_to_int(res[4]);

      // size filter
      if ((end - start) < this->getFilterSize()) {
        continue;
      }

      // roofing the value
      if (value > 5) {
        value = 5;
      }

      // Init the chromosome map 
      if (!(this->data.count(chromosome) > 0)) {
        unordered_map<unsigned int, unordered_map<long, short>> tempMap;
        this->data[chromosome] = tempMap;
        unordered_map<long, short> tempMap2;
        this->data[chromosome][0] = tempMap2;
        this->data[chromosome][1] = tempMap2;
        this->data[chromosome][2] = tempMap2;
        this->data[chromosome][3] = tempMap2;
        this->data[chromosome][4] = tempMap2;
        this->data[chromosome][5] = tempMap2;
      }
      // filling with values
      for (long i = start; i <= end; i++) {
        int value_to_insert = 0; 
        if ((this->data[chromosome][value]).find(i) != (this->data[chromosome][value]).end()) {
          value_to_insert += this->data[chromosome][value][i] + 1; 
        }
        this->data[chromosome][value].insert_or_assign(i, value_to_insert);
      }
    }
    BOOST_LOG_TRIVIAL(info) << " with " << nbLigneFile << " events detected " << endl;
  }
  BOOST_LOG_TRIVIAL(info) << "Ended with " << this->getNbFile() << " files" << endl;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getDataWhole " << endl;
}

/**
 * @brief Method used to commpute counts from data in memory on the whole genome
 * @param none
 * @return none
 **/
void cnvCompare::computeCountsWhole() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::computeCountsWhole " << endl;
  BOOST_LOG_TRIVIAL(info) << "Computing counts" << endl;
  // struct timeval tbegin, tend;
  string ligne;
  string ligneCNV;
  string mot;
  string header = "#";
  int i = 0;
  string currentChr;
  string chromosome;
  string s_type;
  string s_start;
  string s_end;
  string s_value;
  string outFileName;

  map<string, string>::iterator myIterA;
  for (myIterA = this->fileMap.begin(); myIterA != this->fileMap.end(); myIterA++) {
    ligne = myIterA->first;
    string s = myIterA->second;
    if (s == "control") {
      continue;
    }

    ifstream cnvStream(ligne.c_str());

    BOOST_LOG_TRIVIAL(info) << "\tReading file " << ligne << endl;
    fs::path pathObj(ligne);
    if (pathObj.has_extension()) {
      string extension = pathObj.extension().string(); 
      BOOST_LOG_TRIVIAL(trace) << "Extension detected : " << extension << endl;
      BOOST_LOG_TRIVIAL(trace) << "Suffix is : " << this->getSuffix() << endl;
      outFileName = pyReplace(ligne, extension, ("." + this->getSuffix() + extension));
      BOOST_LOG_TRIVIAL(trace) << "outFileName is : " << outFileName << endl;
    } else {
      outFileName = ligne + "." + this->getSuffix();
    }
    if (outFileName == ligne) {
      BOOST_LOG_TRIVIAL(error) << "The output filename " << outFileName << " is the same as the input " << ligne << " : it will replace the original file : Stopping execution" << endl;
      exit(1);
    }

    BOOST_LOG_TRIVIAL(info) << "\t\tWriting in file " << outFileName << endl;
    ofstream outStream;
    outStream.open(outFileName.c_str(), ios::out);
    while (getline(cnvStream, ligneCNV)) {
      // need to deal with header "#"
      if (ligneCNV.find(header) == 0) {
        if (this->getFormat() == "VCF") {
          outStream << ligneCNV << endl;
        }
        continue;
      }
      vector<string> res;
      if (this->getFormat() == "BED") {
        res = this->parseBEDLine(ligneCNV);
      } else {
        res = this->parseVCFLine(ligneCNV);
      }

      // type conversion
      chromosome = res[0];
      s_type = res[3];
      BOOST_LOG_TRIVIAL(trace) << "s_type = " << s_type << endl;
      if ((s_type != "DUP") && (s_type != "DEL")) {
        outStream << ligneCNV << endl;
        continue;
      }

      long start = string_to_int(res[1]);
      long end = string_to_int(res[2]);
      unsigned int value = string_to_int(res[4]);

      if (value > 5) {
        value = 5;
      }

      // counts
      double total = 0;
      for (long i = start; i <= end; i++) {
        total += (double)(this->data[chromosome][value][i]);
      }
      double mean = total / (double)((end-start)+1);

      // need to adapt the output according to the choosen format
      if (this->getFormat() == "BED") {
        outStream << chromosome << "\t" << start << "\t" << end << "\t";
        if (value > 2) {
          outStream << "DUP\t";
        } else {
          outStream << "DEL\t";
        }
        outStream << value << "\t" << mean << "/" << this->getNbFile() << endl;
      } else {
        // output VCF
        res = this->parseVCFLine(ligneCNV);
        istringstream issLigne(ligneCNV);
        istringstream issInfo;
        string mot;
        string info;
        string infomot;
        string svtype;
        string ciend;
        string value;
        i = 0;

        while (getline(issLigne, mot, '\t')) {
          switch (i) {
          case 0:
            outStream << mot;
            break;
          case 7:
            outStream << "\t";
            info = mot;
            issInfo.str(info);
            while (getline(issInfo, infomot, ';')) {
              if (infomot.find("SVTYPE=") == 0) {
                svtype = parseOnSep(infomot, "=")[1];
                continue;
              }
              if (infomot.find("END=") == 0) {
                ciend = parseOnSep(infomot, "=")[1];
                continue;
              }
              if (infomot.find("VALUE=") == 0) {
                value = parseOnSep(infomot, "=")[1];
                continue;
              }
              outStream << infomot << ";";
            }
            outStream << "END=" << ciend << ";VALUE=" << value << ";SVTYPE=";
            if (string_to_int(value) > 2) {
              outStream << "DUP;";
            } else {
              outStream << "DEL;";
            }
            outStream << "COUNT=" << mean << "/" << this->getNbFile();
            break;
          default:
            outStream << "\t" << mot;
            break;
          }
          i++;
        }
        outStream << endl;
      }
    }
    outStream.close();
  }
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::computeCountsWhole " << endl;
}


/**
 * @brief Method used to commpute counts with the fast on the whole genome : default method
 * @param none
 * @return none
 **/
void cnvCompare::computeCountsFast() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::computeCountsFast " << endl;
  BOOST_LOG_TRIVIAL(info) << "Computing counts" << endl;
  // struct timeval tbegin, tend;
  string ligne;
  string ligneCNV;
  string mot;
  string header = "#";
  int i = 0;
  string currentChr;
  string chromosome;
  string s_type;
  string s_start;
  string s_end;
  string s_value;
  string outFileName;

  map<string, string>::iterator myIterA;
  for (myIterA = this->fileMap.begin(); myIterA != this->fileMap.end(); myIterA++) {
    ligne = myIterA->first;
    string s = myIterA->second;
    if (s == "control") {
      continue;
    }
    //managing the output file format
    ifstream cnvStream(ligne.c_str());
    BOOST_LOG_TRIVIAL(info) << "\tReading file " << ligne << endl;
    fs::path pathObj(ligne);
    if (pathObj.has_extension()) {
      string extension = pathObj.extension().string(); 
      BOOST_LOG_TRIVIAL(trace) << "Extension detected : " << extension << endl;
      BOOST_LOG_TRIVIAL(trace) << "Suffix is : " << this->getSuffix() << endl;
      outFileName = pyReplace(ligne, extension, ("." + this->getSuffix() + extension));
      BOOST_LOG_TRIVIAL(debug) << "outFileName is : " << outFileName << endl;
    } else {
      outFileName = ligne + "." + this->getSuffix();
    }
    if (outFileName == ligne) {
      BOOST_LOG_TRIVIAL(error) << "The output filename " << outFileName << " is the same as the input " << ligne << " : it will replace the original file : Stopping execution" << endl;
      exit(1);
    }

    // parsing input file format
    BOOST_LOG_TRIVIAL(info) << "\t\tWriting in file " << outFileName << endl;
    ofstream outStream;
    outStream.open(outFileName.c_str(), ios::out);
    while (getline(cnvStream, ligneCNV)) {
      // need to deal with header "#"
      if (ligneCNV.find(header) == 0) {
        if (this->getFormat() == "VCF") {
          outStream << ligneCNV << endl;
        }
        continue;
      }
      vector<string> res;
      if (this->getFormat() == "BED") {
        res = this->parseBEDLine(ligneCNV);
      } else  {
        res = this->parseVCFLine(ligneCNV);
      }

      // type conversion
      chromosome = res[0];
      s_type = res[3];
      BOOST_LOG_TRIVIAL(debug) << "s_type = " << s_type << endl;
      if ((s_type != "DUP") && (s_type != "DEL")) {
        outStream << ligneCNV << endl;
        continue;
      }

      long start = string_to_int(res[1]);
      long end = string_to_int(res[2]);
      unsigned int value = string_to_int(res[4]);
      // roofing the value
      if (value > 5) {
        value = 5;
      }

      // counts
      BOOST_LOG_TRIVIAL(debug) << "Computing counts " << chromosome << ":" << start << "-" << end << ":" << s_type << value << endl;
      double total = 0;
      map<long, short>::iterator it; 
      short lastValue = 0; 
      long lastPoint = 0; 
      for (it = this->breakpoints[chromosome][value].find(start) ; it != next(this->breakpoints[chromosome][value].find(end), 1) ; ++it) {
        BOOST_LOG_TRIVIAL(debug) << "\toutFileName is : " << outFileName << endl;

        BOOST_LOG_TRIVIAL(debug) << "\tcurrent BP is " << it->first << ":" << it->second;
        if (lastPoint != 0) {
          BOOST_LOG_TRIVIAL(debug) << "\t\tadding " << ((it->first + 1) - lastPoint) * lastValue;
          total += ((it->first + 1) - lastPoint) * lastValue;
          BOOST_LOG_TRIVIAL(debug) << "\t\ttotal is now " << total;
          lastPoint = it->first;
          lastValue = it->second; 
        } else {
          BOOST_LOG_TRIVIAL(debug) << "\t\tNot counting it";
          lastPoint = it->first;
          lastValue = it->second; 
        }
      }
      
      double mean = total / (double)((end-start)+1);
      BOOST_LOG_TRIVIAL(debug) << "\tMean = " << mean; 

      // need to adapt the output according to the choosen format
      if (this->getFormat() == "BED") {
        outStream << chromosome << "\t" << start << "\t" << end << "\t";
        if (value > 2) {
          outStream << "DUP\t";
        } else {
          outStream << "DEL\t";
        }
        outStream << value << "\t" << mean << "/" << this->getNbIndividual() << endl;
      } else {
        // output VCF
        res = this->parseVCFLine(ligneCNV);
        istringstream issLigne(ligneCNV);
        istringstream issInfo;
        string mot;
        string info;
        string infomot;
        string svtype;
        string ciend;
        string value;
        i = 0;

        while (getline(issLigne, mot, '\t')) {
          switch (i) {
          case 0:
            outStream << mot;
            break;
          case 7:
            outStream << "\t";
            info = mot;
            issInfo.str(info);
            while (getline(issInfo, infomot, ';')) {
              if (infomot.find("SVTYPE=") == 0) {
                svtype = parseOnSep(infomot, "=")[1];
                continue;
              }
              if (infomot.find("END=") == 0) {
                ciend = parseOnSep(infomot, "=")[1];
                continue;
              }
              if (infomot.find("VALUE=") == 0) {
                value = parseOnSep(infomot, "=")[1];
                continue;
              }
              outStream << infomot << ";";
            }
            outStream << "END=" << ciend << ";VALUE=" << value << ";SVTYPE=";
            if (string_to_int(value) > 2) {
              outStream << "DUP;";
            } else {
              outStream << "DEL;";
            }
            outStream << "COUNT=" << mean << "/" << this->getNbIndividual();
            break;
          default:
            outStream << "\t" << mot;
            break;
          }
          i++;
        }
        outStream << endl;
      }
    }
    outStream.close();
  }
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::computeCountsFast " << endl;
}


/**
 * @brief Method used to collect data from input files : fast mode
 * @param none
 * @return none
 **/
void cnvCompare::getDataFast() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getDataFast " << endl;
  BOOST_LOG_TRIVIAL(info) << "Gathering data" << endl;
  // struct timeval tbegin, tend;
  string ligne;
  string ligneCNV;
  string mot;
  string header = "#";
  string currentChr;
  string chromosome;
  string s_type;
  string s_start;
  string s_end;
  string s_value;
  vector <short> levelValues(6, 0); 
  

  // tsv parsing
  map<string, string>::iterator myIterA;
  for (myIterA = this->fileMap.begin(); myIterA != this->fileMap.end(); myIterA++) {
    ligne = myIterA->first;
    ifstream cnvStream(ligne.c_str());
    BOOST_LOG_TRIVIAL(info) << "\tReading file " << ligne << "\t";
    this->nbFile++;
    long nbLigneFile = 0;
    short nbOfConcernedIndividual = 0;
    while (getline(cnvStream, ligneCNV)) {
      // need to deal with header "#"
      if (ligneCNV.find(header) == 0) {
        this->watchHeader(ligneCNV);
        continue;
      }
      nbLigneFile++;
      vector<string> res;
      if (this->getFormat() == "BED") {
        res = this->parseBEDLine(ligneCNV);
      } else {
        res = this->parseVCFLine(ligneCNV);

        // check if the vcf parsing was ok
        if (res.size() == 0) {
          BOOST_LOG_TRIVIAL(error) << "\tParsing VCF line : " << ligneCNV << " failed, passing line"; 
          continue; 
        }
        // pass if not del or dup 
        if ((res[3] != "DEL") and (res[3] != "DUP")) {
          BOOST_LOG_TRIVIAL(info) << "\tPassing VCF line : " << ligneCNV << " not DEL nor DUP, passing line"; 
          continue;
        }
        
        nbOfConcernedIndividual = string_to_int(res[5]);
      }

      // type conversion
      chromosome = res[0];
      s_type = res[3];
      long start = string_to_int(res[1]);
      long end = string_to_int(res[2]);
      

      // size filter
      if ((end - start) < this->getFilterSize()) {
        continue;
      }

      // value management
      unsigned int value = string_to_int(res[4]);
      if (value > 5) {
        value = 5;
      }
      
      levelValues[0] = string_to_int(parseOnSep(res[6], ",")[0]);
      levelValues[1] = string_to_int(parseOnSep(res[6], ",")[1]);
      levelValues[2] = string_to_int(parseOnSep(res[6], ",")[2]);
      levelValues[3] = string_to_int(parseOnSep(res[6], ",")[3]);
      levelValues[4] = string_to_int(parseOnSep(res[6], ",")[4]);
      levelValues[5] = string_to_int(parseOnSep(res[6], ",")[5]);
      
      if (value != -1) {
        levelValues[value] = 1;
      }


      for (int cn = 0 ; cn <= 5 ; cn ++) {
        int count = 0; 
        BOOST_LOG_TRIVIAL(debug) << "\tCnv values for this CNV ; cn = " <<  cn  << ", counts = " << levelValues[cn];
        count = levelValues[cn]; 
        if (count == 0) {
          BOOST_LOG_TRIVIAL(debug) << "\t\tnothing to insert";
          continue; 
        }
        
        BOOST_LOG_TRIVIAL(debug) << "\t\twill insert : " << chromosome << ":" << start << "-" << end << " ; cnv : " << cn;
        
        // fill empty map if chr is not existing
        if (!(this->breakpoints.count(chromosome) > 0)) {
          BOOST_LOG_TRIVIAL(debug) << "\t\t\tCreating breakpoint map for this chromosome"; 
          unordered_map<unsigned int, map<long, short> > tempMap;
          this->breakpoints[chromosome] = tempMap;
          map<long, short> tempList;
          this->breakpoints[chromosome][0] = tempList;
          this->breakpoints[chromosome][1] = tempList;
          this->breakpoints[chromosome][2] = tempList;
          this->breakpoints[chromosome][3] = tempList;
          this->breakpoints[chromosome][4] = tempList;
          this->breakpoints[chromosome][5] = tempList;
        }

        // look for the start / end values
        // if the map is empty do not try to browse it, just insert the start and end values and treat the next line. 
        if (this->breakpoints[chromosome][cn].empty()) {
          BOOST_LOG_TRIVIAL(debug) << "\t\t\tMap was empty : so just inserting start & end";
          this->breakpoints[chromosome][cn][start] = 1;
          this->breakpoints[chromosome][cn][end] = 0;
          continue; 
        }

        // insert start and end values in the sorted map
        short lastCount = 0;
        map<long, short>::iterator it, inserted_it, it_beforestart, it_beforeend, it_afterstart, it_afterend;

        // need to get old values
        it_beforestart = breakpoints[chromosome][cn].lower_bound(start);
        it_beforeend = breakpoints[chromosome][cn].lower_bound(end);

        it_afterstart = breakpoints[chromosome][cn].upper_bound(start);
        it_afterend = breakpoints[chromosome][cn].upper_bound(end);
        
        if (it_beforestart != breakpoints[chromosome][cn].end()) {
          it_beforestart --;
          BOOST_LOG_TRIVIAL(debug) << "\t\tBreakpoint before start is " << it_beforestart->first << ":" << it_beforestart->second;
        } else {
          BOOST_LOG_TRIVIAL(debug) << "\t\tBreakpoint before start is after the current end of the map";
        }
        if (it_beforeend != breakpoints[chromosome][cn].end()) {
          it_beforeend --;
          BOOST_LOG_TRIVIAL(debug) << "\t\tBreakpoint before end is " << it_beforeend->first << ":" << it_beforeend->second;
        } else {
          BOOST_LOG_TRIVIAL(debug) << "\t\tBreakpoint before end is after the current end of the map";
        }
        if (it_afterstart != breakpoints[chromosome][cn].end()) {
          BOOST_LOG_TRIVIAL(debug) << "\t\tBreakpoint after start is " << it_afterstart->first << ":" << it_afterstart->second;
        } else {
          BOOST_LOG_TRIVIAL(debug) << "\t\tBreakpoint after start is after the current end of the map";
        }
        if (it_afterend != breakpoints[chromosome][cn].end()) {
          BOOST_LOG_TRIVIAL(debug) << "\t\tBreakpoint after end is " << it_afterend->first << ":" << it_afterend->second;
        } else {
          BOOST_LOG_TRIVIAL(debug) << "\t\tBreakpoint after end is after the current end of the map";
        }

        // manage begin of the map
        if ((it_beforestart == breakpoints[chromosome][cn].begin()) || (it_beforestart == breakpoints[chromosome][cn].end())){
          lastCount = 0;
        } else {
          lastCount =  it_beforestart->second;
        }

        // insert start point 
        breakpoints[chromosome][cn].insert_or_assign(start, lastCount + count);
        BOOST_LOG_TRIVIAL(debug) << "\t\tStart inserted " << start << ":" << lastCount + count << " at cn level " << cn;

        // get last value of the interval & manage end of the map
        if (it_beforeend == breakpoints[chromosome][cn].end()) {
          lastCount = 0;
        } else {
          lastCount =  it_beforeend->second;
        }

        // modify all value until end
        for (it = it_afterstart ; it != it_afterend ; ++ it) {
          if (it != breakpoints[chromosome][cn].end()) {
            breakpoints[chromosome][cn][it->first] += count;
            BOOST_LOG_TRIVIAL(debug) << "\t\tChanging breakpoints " << it->first << ":" << breakpoints[chromosome][cn][it->first] - count << " to " << breakpoints[chromosome][cn][it->first] << " at cn level " << cn;
          }
        }

        // insert the end 
        breakpoints[chromosome][value].insert_or_assign(end, lastCount);
        BOOST_LOG_TRIVIAL(debug) << "\t\tEnd inserted " << end << ":" << lastCount << " at cn level " << cn;
        BOOST_LOG_TRIVIAL(debug) << "\t\tSize of breakpoints at chr : " << chromosome << " and value " << value << " : " << breakpoints[chromosome][cn].size();

        // a count for large files to be sure that everything went well
        if ((nbLigneFile % 10000) == 0) {
          BOOST_LOG_TRIVIAL(info) << "\t" << nbLigneFile << " events detected, still in progress" << endl;
        }
      }
    }
    BOOST_LOG_TRIVIAL(info) << " with " << nbLigneFile << " events detected " << endl;
  }
  BOOST_LOG_TRIVIAL(info) << "Ended with " << this->getNbFile() << " files" << endl;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getDataFast " << endl;
}


/**
 * @brief Getter for the number of input files
 * @param none
 * @return integer value : the number of file from the instance of the cnvCompare class
 **/
short int cnvCompare::getNbFile() { 
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getNbFile " << endl;
  return this->nbFile; 
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getNbFile " << endl;
}


/**
 * @brief Getter for the name of the input file list
 * @param none
 * @return string value : the name of the input file list
 **/
string cnvCompare::getInputFile() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getInputFile " << endl;
  return this->inputFile;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getInputFile " << endl;
}


/**
 * @brief Getter for the name of the control file list
 * @param none
 * @return string value : the name of the control file list
 **/
string cnvCompare::getControlFile() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getControlFile " << endl;
  return this->controlFile;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getControlFile " << endl;
}


/**
 * @brief Method used to fill the file map
 * @param incFile : string value : the path to the file containing a list of  files 
 * @param status : string value : precising if the files are controls or inputs
 * @return int value : a return code : 0 if ok, > 0 if not ok
 **/
int cnvCompare::fillMap(string incFile, string status) {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::fillMap " << endl;
  ifstream allStream(incFile.c_str());
  string ligne;
  if (allStream) {
    while (getline(allStream, ligne)) {
      ifstream cnvStream(ligne.c_str());
      BOOST_LOG_TRIVIAL(info) << "\tadding file " << ligne << " as " << status << endl;
      if (!IsFileReadable(ligne)) {
        BOOST_LOG_TRIVIAL(warning) << "File provided as " << status << " : " << ligne << " is not accessible : passsing file" << endl;
        continue;
      }
      this->nbFile++;
      this->fileMap[ligne] = status;
    }
  }
  
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::fillMap " << endl;
  return this->nbFile;
}


/**
 * @brief Getter for filter size value
 * @param none
 * @return int value : the filter size value
 **/
int cnvCompare::getFilterSize() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getFilterSize " << endl;
  return this->filterSize;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getFilterSize " << endl;
}


/**
 * @brief Setter for the file format
 * @param incVCFChoice A boolean indicating if the files are in VCF format
 * @param incBEDChoice A boolean indicating if the files are in BED format
 * @return none
 **/
void cnvCompare::setFormat(bool incVCFChoice, bool incBEDChoice) {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::setFormat " << endl;
  this->useVCFFormat = incVCFChoice;
  this->useBEDFormat = incBEDChoice;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::setFormat " << endl;
}


/**
 * @brief Getter for the file formate
 * @param none
 * @return string value : the file format
 **/
string cnvCompare::getFormat() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getFormat " << endl;
  if (this->useVCFFormat) {
    return "VCF";
  } else {
    return "BED";
  }
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getFormat " << endl;
}


/**
 * @brief Method used to populate the chromosome map
 * @param none
 * @return int value : a return code : 0 if ok, > 0 if not ok
 **/
int cnvCompare::populateChr() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::populateChr " << endl;
  this->chromosomeMap.insert(this->chromosomeMap.end(), {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"});
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::populateChr " << endl;
  return 0; 
}


/**
 * @brief Setter used to define the suffix to use
 * @param incSuffix : string defining the suffix value
 * @return none
 **/
void cnvCompare::setSuffix(string incSuffix) {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::setSuffix " << endl;
  this->suffix = incSuffix;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::setSuffix " << endl;
}


/**
 * @brief getter used to get the suffix to use
 * @param none
 * @return string value : defining the suffix value
 **/
string cnvCompare::getSuffix() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getSuffix " << endl;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getSuffix " << endl;
  return this->suffix; 
}


/**
 * @brief Setter used to define the dict file to use to populate the chromosome map
 * @param incDictFile : string defining the dict file path
 * @return none
 **/
void cnvCompare::setDictFile(string incDictFile) {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::setDictFile " << endl;
  this->dictFile = incDictFile;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::setDictFile " << endl;
}


/**
 * @brief getter used to get the dict file to use
 * @param none
 * @return string value defining the dict file
 **/
string cnvCompare::getDictFile() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getDictFile " << endl;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getDictFile " << endl;
  return this->dictFile; 
}


/**
 * @brief Method used to parse a dictionnary file and to fill the apporpriate container. 
 * @param incDictFile : string defining path to a dictionnary file
 * @return none
 **/
void cnvCompare::parseDictFile(string incDictFile) {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::parseDictFile " << endl;
  string ligne;
  string mot;
  string header = "@SQ";

  // Parsing dict file
  BOOST_LOG_TRIVIAL(info) << "Reading Dict file " << incDictFile << endl;
  ifstream dictStream(incDictFile.c_str());
  while (getline(dictStream, ligne)) {
    string chrName; 
    // need to deal with header 
    if (ligne.find(header) == 0) {
      vector<string> vline = parseOnSep("\t", ligne); 
      vector<string>::iterator myIter; 
      for (myIter = vline.begin(); myIter != vline.end(); myIter ++) {
        if ((*myIter).find("SN") == 0) {
          chrName = parseOnSep(":", (*myIter))[1]; 
        }
      }
      BOOST_LOG_TRIVIAL(info) << "adding " << chrName << " to the chromosome map" << endl;
      this->chromosomeMap.push_back(chrName);
    }
  }
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::parseDictFile " << endl;
}

/**
 * @brief Setter used to define if a dict file is present
 * @param incBool : boolean indicating if a dict file is present
 * @return none
 **/
void cnvCompare::setHasDict(bool incBool) {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::setHasDict " << endl;
  this->hasDict = incBool;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::setHasDict " << endl;
}


/**
 * @brief Getter for number of individuals (different from number of files)
 * @param none
 * @return short value : the number of individual in memory
 **/
short int cnvCompare::getNbIndividual() {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getNbIndividual " << endl;
  return this->nbIndividual;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getNbIndividual " << endl;
}

/**
 * @brief Method used to get the number of individuals on a header line for a VCF
 * @param incLine : A string containing the header line
 * @return none
 **/
void cnvCompare::watchHeader(string incLine) {
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::watchHeader " << endl;
  string goodHeader = "#CHR"; 
  short initIndiv = this->getNbIndividual();
  if (incLine.find(goodHeader) == 0) {
    vector<string> lineTable = parseOnSep(incLine, "\t");
    for (long unsigned int n = 9 ; n < lineTable.size() ; n ++) {
      this->nbIndividual += 1;
    }
    BOOST_LOG_TRIVIAL(info) << "Added " << this->getNbIndividual() - initIndiv << " individuals to the list, total is now : " << this->getNbIndividual() << endl;
  }

  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::watchHeader" << endl;
}
