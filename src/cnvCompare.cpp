// C++ std libs
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <map>
#include <utility>
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
#include <cmath>
//#include <ranges>

// Boost 
#include <boost/filesystem/path.hpp>

// plog 
#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>


// utils
#include "utils.h"

// Classes
#include "cnvCompare.h"

// namespaces
using namespace std;
using namespace boost;
namespace fs = boost::filesystem;

/**
 * @brief default constructor (useless)
 * @param none
 * @return none
 **/
cnvCompare::cnvCompare() {
  PLOG(plog::verbose) << "Entering cnvCompare::cnvCompare default constructor";
  PLOG(plog::verbose) << "Leaving cnvCompare::cnvCompare default constructor";
}


/**
 * @brief Constructor of the class without controls
 * @param iF file containing the list of input files
 * @param iF number of threads 
 * @param s size of the filter to apply
 * @return none
 **/
cnvCompare::cnvCompare(string iF, int nT, int s) {
  PLOG(plog::verbose) << "Entering cnvCompare::cnvCompare (string, int, int) Ctor";
  PLOG(plog::verbose) <<  "cnvCompare constructor called";
  PLOG(plog::verbose) <<  "sizeof(char) = " << sizeof(char);
  PLOG(plog::verbose) <<  "sizeof(string) = " << sizeof(string);
  PLOG(plog::verbose) <<  "sizeof(int32_t) = " << sizeof(int32_t);
  PLOG(plog::verbose) <<  "sizeof(int) = " << sizeof(int);
  
  // init of the needed class members
  this->inputFile = iF;
  this->nbThread = nT;
  this->useControls = false;
  this->filterSize = s;
  this->suffix = "count";
  // Adding sample files 
  int n = this->fillMap(iF, "sample");
  PLOG(plog::info) << n << " files added to sample list";
  // populating chromosome either by parsing the dict of with a default human map
  if (this->hasDict) {
    this->parseDictFile(this->getDictFile());
    
  } else {
    int ret = this->populateChr();
    if (ret != 0) {
      PLOG(plog::error) << "Chromosome population couldn't be filled in";
    }
  }
  PLOG(plog::verbose) << "Leaving cnvCompare::cnvCompare (string, int, int) Ctor"; 
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
  PLOG(plog::verbose) << "Entering cnvCompare::cnvCompare (string, string, int, int) Ctor";
  // init of the needed class members
  this->inputFile = iF;
  this->controlFile = cF;
  this->nbThread = nT;
  this->useControls = true;
  this->filterSize = s;
  this->suffix = "count";
  // Adding sample files & control files
  int n = this->fillMap(iF, "sample");
  PLOG(plog::info) << n << " files added to sample list";
  int m = this->fillMap(cF, "control");
  PLOG(plog::info) << m << " files added to control list";
  // populating chromosome either by parsing the dict of with a default human map
  if (this->hasDict) {
    this->parseDictFile(this->getDictFile());
  } else {
    int ret = this->populateChr();
    if (ret != 0) {
      PLOG(plog::error) << "Chromosome population couldn't be filled in";
    }
  }
  PLOG(plog::verbose) << "Leaving cnvCompare::cnvCompare (string, string, int, int) Ctor";
}


/**
 * @brief Main loop used to get data, and compute counts chr by chr
 * @param none
 * @return none
 **/
void cnvCompare::mainLoop() {
  PLOG(plog::verbose) << "Entering cnvCompare::mainLoop ";
  for (auto &a : this->chromosomeMap) {
    // getting data
    this->getDatabyChr(a);
    // computing counts 
    this->computeCountsbyChr(a);
    // GC 
    this->cleanData();
  }
  PLOG(plog::verbose) << "Leaving cnvCompare::mainLoop ";
}

/**
 * @brief Alt loop used to get data, and compute counts on the whole genome ; needs huge amount of RAM
 *        Warning : this method is deprecated. Do not use. 
 * @param none
 * @return none
 **/
void cnvCompare::altLoop() {
  PLOG(plog::verbose) << "Entering cnvCompare::altLoop ";
  // getting data
  this->getDataWhole();
  // computing counts 
  this->computeCountsWhole();
  // GC 
  this->cleanData();
  PLOG(plog::verbose) << "Leaving cnvCompare::altLoop ";
}

/**
 * @brief fast loop used to get data, and compute counts on intervals ; experimental feature
 * @param none
 * @return none
 **/
void cnvCompare::fastLoop() {
  PLOG(plog::verbose) << "Entering cnvCompare::fastLoop ";
  // getting data
  this->getDataFast();
  // computing counts 
  this->computeCountsFast();
  // GC 
  this->cleanData();
  PLOG(plog::verbose) << "Leaving cnvCompare::fastLoop ";
}


/**
 * @brief Method used to emptying data for the chromosomes
 * @param none
 * @return none
 **/
void cnvCompare::cleanData() { 
  PLOG(plog::verbose) << "Entering cnvCompare::cleanData ";
  this->dataByChr.clear();
  this->data.clear();
  this->breakpoints.clear(); 
  PLOG(plog::verbose) << "Leaving cnvCompare::cleanData "; 
}

/**
 * @brief Method used to collect data from input files
 * @param incChr the chr to filter
 * @return none
 **/
void cnvCompare::getDatabyChr(string incChr) {
  PLOG(plog::verbose) << "Entering cnvCompare::getDatabyChr ";
  PLOG(plog::info) << "Gathering data for chr " << incChr;
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
    PLOG(plog::info) << "\tReading file " << ligne << "\t";
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
    PLOG(plog::info) << nbLigneFile << " events detected ";
  }
  PLOG(plog::info) << "Ended with " << this->getNbFile() << " files";
  PLOG(plog::verbose) << "Leaving cnvCompare::getDatabyChr ";
}


/**
 * @brief Method used to compute counts from in memory data chr by chr
 * @param incChr the chr to filter
 * @return none
 **/
void cnvCompare::computeCountsbyChr(string incChr) {
  PLOG(plog::verbose) << "Entering cnvCompare::computeCountsbyChr ";
  std::cout.precision(3);
  PLOG(plog::info) << "Computing counts for chr " << incChr;
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
    PLOG(plog::info) << "\tReading file " << ligne;

    // managing output file names
    fs::path pathObj(ligne);
    if (pathObj.has_extension()) {
      string extension = pathObj.extension().string(); 
      PLOG(plog::debug) << "Extension detected : " << extension;
      outFileName = pyReplace(ligne, extension, "." + this->getSuffix() + extension);
    } else {
      outFileName = ligne + "." + this->getSuffix();
    }
    if (outFileName == ligne) {
      PLOG(plog::error) << "The output filename " << outFileName << " is the same as the input " << ligne << " : it will replace the original file : Stopping execution";
      exit(1);
    }
    PLOG(plog::info) << "\t\tWriting in file " << outFileName << " in format : " << this->getFormat();
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
      PLOG(plog::verbose) << "s_type = " << s_type;
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
  PLOG(plog::verbose) << "Leaving cnvCompare::computeCountsbyChr ";
}


/**
 * @brief Method used to parse a line from a BED file
 * @param incLine the line to parse
 * @return vector<string> containing values : chr, start, end, type, value
 * @todo need to check that all the values are present
 **/
// output a vector containing : chr start end type value from a BED line
vector<string> cnvCompare::parseBEDLine(string incLine) {
  PLOG(plog::verbose) << "Entering cnvCompare::parseBEDLine ";
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
  PLOG(plog::verbose) << "Leaving cnvCompare::parseBEDLine ";
}


/**
 * @brief Method used to parse a line from a VCF file
 * @param incLine the line to parse
 * @return vector<string> containing values : chr, start, end, type, value
 * @todo need to check that all the values are present
 **/
vector<string> cnvCompare::parseVCFLine(string incLine) {
  PLOG(plog::verbose) << "Entering cnvCompare::parseVCFLine ";
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
  vector<short> counts(7, 0); 
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
        if ((infomot.find("VALUE=") == 0) or (infomot.find("CN=") == 0)){
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
      for (long unsigned int n = 0 ; n < GTInfo.size() ; n ++) {
        if (strcmp(GTInfo[n].c_str(), "GT") == 0) {
          GTindex = n;
          PLOG(plog::debug) << "GT index was found : " << GTindex; 
        }
        if (! valueFound) {
          if (strcmp(GTInfo[n].c_str(), "CN") == 0) {
            CNindex = n;
            PLOG(plog::debug) << "CN index was found : " << GTindex; 
          }
        }
      }
      ++i;
      continue;
    }
    
      // getting the non wild indiv 
    if (i >= 9) {
      // first checking if the GT field has been found
      if (GTindex == -1) {
        i++; 
        continue;
      }

      if ((temp["SVTYPE"] != "DEL") and (temp["SVTYPE"] != "DUP") and (temp["SVTYPE"] != "INV") and (temp["SVTYPE"] != "CNV")){
        break; 
      }

      // copy number value 
      PLOG(plog::debug) << "\tComputing cn value";  
      if ((! valueFound) && (CNindex == -1)) {
        PLOG(plog::info) << "\t\tNo copy number value found on the VCF line " << incLine;
        PLOG(plog::info) << "\t\tPlease check the VCF specifications. Will try to infer it with GT field";
        if (GTindex != -1) {
          if ((temp["SVTYPE"] != "INV") && (temp["SVTYPE"] != "CNV")) {
            string GT = parseOnSep(mot, ":")[GTindex];
            PLOG(plog::debug) << "\t\tGT is " << GT;
            CNValue_i = this->inferCNfromGT(GT, temp["SVTYPE"]);
            if (CNValue_i > 5) {
              CNValue_i = 5;
            }
          } 
        } else {
            PLOG(plog::info) << "\t\tNo GT field : passing";
            passGT = true;
            break;
        }
        
      } else {
        CNValue_i = string_to_int(temp["VALUE"]);
      }
      PLOG(plog::debug) << "\tCn value is " << CNValue_i;  

 
      // counts number of individual
      PLOG(plog::debug) << "\tCounting individuals";  
      if (GTindex == -1) {
        PLOG(plog::info) << "No GT Found on the VCF line : counting only 1";
        nbOfConcernedIndiv = 1;
        passGT = true;
        break;
      } else {
        PLOG(plog::info) << "GT Found on the VCF line : counting regarding GT";
        string GT = parseOnSep(mot, ":")[GTindex];
        PLOG(plog::info) << "GT is " << GT;
        if ((temp["SVTYPE"] != "INV") && (temp["SVTYPE"] != "CNV")) {   
          CNValue_i = this->inferCNfromGT(GT, temp["SVTYPE"]);
          // not interested in cnv at n=2
          if (CNValue_i == 2) {
            PLOG(plog::info) << "CN is 2 : passing";
            ++i; 
            passGT = true;
            continue; 
          }
          if (CNValue_i > 5) {
            CNValue_i = 5;
          }
          counts[CNValue_i] += 1;
        } else { // inversion
          if (GT != "0/0") {
            counts[6] += 1;
          }
        }
        nbOfConcernedIndiv += 1;
        PLOG(plog::debug) << "\t\tadding 1 concerned individual";;
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
  PLOG(plog::debug) << "\tNumber of concerned individual is " << nbOfConcernedIndiv;
  // transforming the counts vector into string 
  output.push_back(int_to_string(counts[0]) + "," + int_to_string(counts[1]) + "," + int_to_string(counts[2]) + "," + int_to_string(counts[3]) + "," + int_to_string(counts[4]) + "," + int_to_string(counts[5]) + "," + int_to_string(counts[6]));


  PLOG(plog::debug) << "\tWill return output from VCF line : ";
  vector <string>::iterator myIter; 
  for (myIter = output.begin() ; myIter != output.end() ; myIter++ ) {
    PLOG(plog::debug) << "\t\t" << *myIter << endl; 
  }
  PLOG(plog::debug) << "\n";

  PLOG(plog::verbose) << "Leaving cnvCompare::parseVCFLine ";
  return output;
}

/**
 * @brief Method used to parse a line from a VCF file for TRN usage only
 * @param incLine the line to parse
 * @return None
 * @todo need to check that all the values are present
 **/
void cnvCompare::parseVCFLineTRN(string incLine) {
  PLOG(plog::verbose) << "Entering cnvCompare::parseVCFLineTRN ";
  map<string, string> temp;
  string mot;
  short int i = 0;
  temp["SVTYPE"] = "TRN";
  unsigned int trnAPosition = 0;
  unsigned int trnBPosition = 0;
  string trnAChr = ""; 
  string trnBChr = "";
  string trnBulk = ""; 

  // get information from the line
  istringstream issLigne(incLine);
  i = 0;
  while (getline(issLigne, mot, '\t')) {
    if (i == 0) {
      trnAChr = mot;
      ++i;
      continue;
    }
    if (i == 1) {
      trnAPosition = string_to_int(mot);
      ++i;
      continue;
    }

    if (i == 4) {
      trnBulk = mot;
      ++i;
      continue;
    }

    if (i > 4) {
      break;
    }

    ++i;
  }
  PLOG(plog::verbose) << "\ttrnAChr : " << trnAChr;
  PLOG(plog::verbose) << "\ttrnAPosition : " << trnAPosition;
  PLOG(plog::verbose) << "\ttrnBulk : " << trnBulk;

  // modify the second breakpoint if [ or ] detected 
  if (trnBulk.find("[") != string::npos) {
    PLOG(plog::debug) << "\tFound a [ in the BND second breakpoint : " << trnBulk << " at pos " << trnBulk.find("[");
    vector<string> tempBS;
    tempBS =  parseOnSep(trnBulk, "[");
    if (tempBS.size() > 1) {
      trnBulk = tempBS[1];
    }
    PLOG(plog::debug) << "\t\tTurned into : " << trnBulk;
  }
  if (trnBulk.find("]") != string::npos) {
    PLOG(plog::debug) << "\tFound a ] in the BND second breakpoint : " << trnBulk << " at pos : " << trnBulk.find("]");
    vector<string> tempBS;
    tempBS =  parseOnSep(trnBulk, "]");
    if (tempBS.size() > 1) {
      trnBulk = tempBS[1];
    }
    PLOG(plog::debug) << "\t\tTurned into : " << trnBulk;
  }



 
  // get the end point 
  vector<string> trnBulkInfo;
  trnBulkInfo = parseOnSep(trnBulk, ":");
  trnBChr = trnBulkInfo[0];
  trnBPosition = string_to_int(trnBulkInfo[1]);

  // check if there is already a breakpoint existing with these coordinates and get index, insert it if needed
  PLOG(plog::debug) << "\tGetting index A :";
  int indexA = this->insertTrnBreakpoints(trnAChr, trnAPosition);
  PLOG(plog::debug) << "\t\t" << indexA;
  PLOG(plog::debug) << "\tGetting index B :";
  int indexB = this->insertTrnBreakpoints(trnBChr, trnBPosition);
  PLOG(plog::debug) << "\t\t" << indexB;

  // swap index if A > B
  int tmpIndex; 
  if (indexA > indexB) {
    tmpIndex = indexB;
    indexB = indexA;
    indexA = tmpIndex;
  }

  // check if association exists 
  string k = int_to_string(indexA) + "_" + int_to_string(indexB);
  auto it = this->trnAssociation.find(k);
  if (it == this->trnAssociation.end()) {
    PLOG(plog::debug) << "Association : Index A = " << indexA << " : Index B = " << indexB << " was not found, inserting association "; 
    this->trnAssociation[k] = 1;
  } else {
    PLOG(plog::debug) << "Association : Index A = " << indexA << " : Index B = " << indexB << " was found, incrementing "; 
    this->trnAssociation[k] += 1;
  }



  // if (!(this->trnAssociation.count(indexA) > 0)) {
  //   PLOG(plog::debug) << "Index A : " << indexA << " was not found, inserting association "; 
  //   map<int, int> tmpMap;
  //   tmpMap[indexB] = 1;
  //   this->trnAssociation[indexA] = tmpMap;
  //   PLOG(plog::debug) << "Inserting new TRN association : " << indexA << " : " << indexB << ", count = " << 1;
  // } else {
  //   PLOG(plog::debug) << "Index A : " << indexA << " was found"; 
  //   if (!(this->trnAssociation[indexA].count(indexB) > 0)) {
  //     PLOG(plog::debug) << "\tindexB was not found : " << indexB;
  //     map<int, int> tmpMap;
  //     tmpMap[indexB] = 1;
  //     this->trnAssociation[indexA] = tmpMap;
  //     PLOG(plog::debug) << "Inserting new TRN association : " << indexA << " : " << indexB << ", count = " << 1;
  //   } else {
  //     this->trnAssociation[indexA][indexB] ++;
  //     PLOG(plog::debug) << "\tIndexB was found : " << indexB << ", incrementing counts  : " << this->trnAssociation[indexA][indexB];
  //   }
  // }
  PLOG(plog::verbose) << "Leaving cnvCompare::parseVCFLineTRN ";
}

/**
 * @brief Method used to get number of association of TRN in a line from a VCF file
 * @param incLine the line to parse
 * @return the number of association
 * @todo 
 **/
int cnvCompare::getTRNAssoc(string incLine) {
  PLOG(plog::verbose) << "Entering cnvCompare::getTRNAssoc ";
  map<string, string> temp;
  string mot;
  short int i = 0;
  temp["SVTYPE"] = "TRN";
  unsigned int trnAPosition = 0;
  unsigned int trnBPosition = 0;
  string trnAChr = ""; 
  string trnBChr = "";
  string trnBulk = ""; 

  // get information from the line
  istringstream issLigne(incLine);
  i = 0;
  while (getline(issLigne, mot, '\t')) {
    if (i == 0) {
      trnAChr = mot;
      ++ i;
      continue;
    }
    if (i == 1) {
      trnAPosition = string_to_int(mot);
      ++i;
      continue;
    }

    if (i == 4) {
      trnBulk = mot;
      ++i;
      continue;
    }
    if (i > 4) {
      break;
    }
    ++i; 
  }

  // modify the second breakpoint if [ or ] detected 
  if (trnBulk.find("[") != string::npos) {
    PLOG(plog::debug) << "\tFound a [ in the BND second breakpoint : " << trnBulk << " at pos " << trnBulk.find("[");
    vector<string> tempBS;
    tempBS =  parseOnSep(trnBulk, "[");
    if (tempBS.size() > 1) {
      trnBulk = tempBS[1];
    }
    PLOG(plog::debug) << "\t\tTurned into : " << trnBulk;
  }
  if (trnBulk.find("]") != string::npos) {
    PLOG(plog::debug) << "\tFound a ] in the BND second breakpoint : " << trnBulk << " at pos : " << trnBulk.find("]");
    vector<string> tempBS;
    tempBS =  parseOnSep(trnBulk, "]");
    if (tempBS.size() > 1) {
      trnBulk = tempBS[1];
    }
    PLOG(plog::debug) << "\t\tTurned into : " << trnBulk;
  }

  // get the end point 
  vector<string> trnBulkInfo;
  trnBulkInfo = parseOnSep(trnBulk, ":");
  trnBChr = trnBulkInfo[0];
  trnBPosition = string_to_int(trnBulkInfo[1]);

  // check if there is already a breakpoint existing with these coordinates and get index, insert it if needed 
  PLOG(plog::debug) << "\tLooking for : " << trnAChr << ":" << trnAPosition;
  int indexA = this->getTRNIndex(trnAChr, trnAPosition);
  PLOG(plog::debug) << "\tLooking for : " << trnBChr << ":" << trnBPosition;
  int indexB = this->getTRNIndex(trnBChr, trnBPosition);

  // swap index if A > B
  int tmpIndex; 
  if (indexA > indexB) {
    tmpIndex = indexB;
    indexB = indexA;
    indexA = tmpIndex;
  }
  string k = int_to_string(indexA) + "_" + int_to_string(indexB);
  short res = this->trnAssociation[k];
  if (res < 1) {
    res = 1;
  }
  PLOG(plog::verbose) << "Leaving cnvCompare::getTRNAssoc ";
  return res;
}

/**
 * @brief Method used to collect data from input files : Whole mode
 * @param none
 * @return none
 **/
void cnvCompare::getDataWhole() {
  PLOG(plog::verbose) << "Entering cnvCompare::getDataWhole ";
  PLOG(plog::info) << "Gathering data";
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
    PLOG(plog::info) << "\tReading file " << ligne << "\t";
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
    PLOG(plog::info) << " with " << nbLigneFile << " events detected ";
  }
  PLOG(plog::info) << "Ended with " << this->getNbFile() << " files";
  PLOG(plog::verbose) << "Leaving cnvCompare::getDataWhole ";
}

/**
 * @brief Method used to commpute counts from data in memory on the whole genome
 * @param none
 * @return none
 **/
void cnvCompare::computeCountsWhole() {
  PLOG(plog::verbose) << "Entering cnvCompare::computeCountsWhole ";
  PLOG(plog::info) << "Computing counts";
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

    PLOG(plog::info) << "\tReading file " << ligne;
    fs::path pathObj(ligne);
    if (pathObj.has_extension()) {
      string extension = pathObj.extension().string(); 
      PLOG(plog::verbose) << "Extension detected : " << extension;
      PLOG(plog::verbose) << "Suffix is : " << this->getSuffix();
      outFileName = pyReplace(ligne, extension, ("." + this->getSuffix() + extension));
      PLOG(plog::verbose) << "outFileName is : " << outFileName;
    } else {
      outFileName = ligne + "." + this->getSuffix();
    }
    if (outFileName == ligne) {
      PLOG(plog::error) << "The output filename " << outFileName << " is the same as the input " << ligne << " : it will replace the original file : Stopping execution" << endl;
      exit(1);
    }

    PLOG(plog::info) << "\t\tWriting in file " << outFileName;
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
      PLOG(plog::verbose) << "s_type = " << s_type;
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
        outStream << value << "\t" << floor(mean) << "/" << this->getNbFile() << endl;
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
            outStream << "COUNT=" << floor(mean) << "/" << this->getNbFile();
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
  PLOG(plog::verbose) << "Leaving cnvCompare::computeCountsWhole ";
}


/**
 * @brief Method used to commpute counts with the fast on the whole genome : default method
 * @param none
 * @return none
 **/
void cnvCompare::computeCountsFast() {
  PLOG(plog::verbose) << "Entering cnvCompare::computeCountsFast ";
  PLOG(plog::info) << "Computing counts";
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
    PLOG(plog::info) << "\tReading file " << ligne;
    fs::path pathObj(ligne);
    if (pathObj.has_extension()) {
      string extension = pathObj.extension().string(); 
      PLOG(plog::verbose) << "Extension detected : " << extension;
      PLOG(plog::verbose) << "Suffix is : " << this->getSuffix();
      outFileName = pyReplace(ligne, extension, ("." + this->getSuffix() + extension));
      PLOG(plog::debug) << "outFileName is : " << outFileName;
    } else {
      outFileName = ligne + "." + this->getSuffix();
    }
    if (outFileName == ligne) {
      PLOG(plog::error) << "The output filename " << outFileName << " is the same as the input " << ligne << " : it will replace the original file : Stopping execution";
      exit(1);
    }

    // parsing input file format
    PLOG(plog::info) << "\t\tWriting in file " << outFileName;
    ofstream outStream;
    outStream.open(outFileName.c_str(), ios::out);
    while (getline(cnvStream, ligneCNV)) {
      // need todeal with header "#"
      PLOG(plog::debug) << "### NEW LINE ###";
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
      PLOG(plog::debug) << "s_type = " << s_type;
      if ((s_type == "TRN") || (s_type == "BND"))  {
        int nbAssoc;
        nbAssoc = this->getTRNAssoc(ligneCNV);
        // write 
        if (this->getFormat() == "BED") {
          outStream << ligneCNV << endl;
          continue;
        } else {
          // output VCF
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
              outStream << "END=" << ciend << ";VALUE=" << value << ";SVTYPE=BND";
              outStream << ";COUNT=" << nbAssoc << "/" << this->getNbIndividual();
              break;
            default:
              outStream << "\t" << mot;
              break;
            }
            i++;
          }
            outStream << endl;
        }
        continue;
      }


      if ((s_type != "DUP") && (s_type != "DEL") && (s_type != "INV") && (s_type != "CNV")) {
        PLOG(plog::info) << "Found an event type not managed : " << s_type;
        outStream << ligneCNV << endl;
        continue;
      }
      unsigned int value;
      long start = string_to_int(res[1]);
      long end = string_to_int(res[2]);
      if (res[4] == "-1") {
        value = 2;
      } else { 
        value = string_to_int(res[4]);
      }
      // roofing the value
      if (value > 5) {
        value = 5;
      }
      if (s_type == "INV") {
        value = 6;
      }

      // counts
      PLOG(plog::debug) << "Computing counts " << chromosome << ":" << start << "-" << end << ":" << s_type << value;
      double total = 0;
      map<long, short>::iterator it; 
      map<long, short>::iterator endIt;
      short lastValue = 0; 
      long lastPoint = 0; 
      PLOG(plog::debug) << "\toutFileName is : " << outFileName;
      it = this->breakpoints[chromosome][value].find(start);
      PLOG(plog::debug) << "\titerator pointing first to : " << it->first << ":" << it->second;
      it = this->breakpoints[chromosome][value].find(end);
      PLOG(plog::debug) << "\titerator pointing end to : " << it->first << ":" << it->second;

      // manage event on the whole chromosome 
      if (it == this->breakpoints[chromosome][value].end()) {
        endIt = it;
      } else {
        endIt = next(this->breakpoints[chromosome][value].find(end), 1);
      }


      for (it = this->breakpoints[chromosome][value].find(start) ; it != endIt ; ++it) {
        PLOG(plog::debug) << "\tcurrent BP is " << it->first << ":" << it->second;
        if (lastPoint != 0) {
          PLOG(plog::debug) << "\t\tadding " << ((it->first + 1) - lastPoint) * lastValue;
          total += ((it->first + 1) - lastPoint) * lastValue;
          PLOG(plog::debug) << "\t\ttotal is now " << total;
          lastPoint = it->first;
          lastValue = it->second; 
        } else {
          PLOG(plog::debug) << "\t\tNot counting it";
          lastPoint = it->first;
          lastValue = it->second; 
        }
      }
      
      double mean;
      if (total == 0) {
        mean = 1.0;
      } else {      
        mean = total / (double)((end-start)+1);
      }
      PLOG(plog::debug) << "\tMean = " << mean; 

      // need to adapt the output according to the choosen format
      if (this->getFormat() == "BED") {
        outStream << chromosome << "\t" << start << "\t" << end << "\t";
        if (value == 6) {
           outStream << "INV\t";
        } else { 
          if (value > 2) {
          outStream << "DUP\t";
          } else {
            outStream << "DEL\t";
          }
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
            if (value == ".") {
              value = "6";
            }
            if (string_to_int(value) == 6) {
              outStream << "INV;";
            } else { 
              if (string_to_int(value) > 2) {
                  outStream << "DUP;";
                } else {
                  outStream << "DEL;";
                }
            }
            outStream << "COUNT=" << floor(mean) << "/" << this->getNbIndividual();
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
  PLOG(plog::verbose) << "Leaving cnvCompare::computeCountsFast ";
}


/**
 * @brief Method used to collect data from input files : fast mode
 * @param none
 * @return none
 **/
void cnvCompare::getDataFast() {
  PLOG(plog::verbose) << "Entering cnvCompare::getDataFast ";
  PLOG(plog::info) << "Gathering data";
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
    PLOG(plog::info) << "\tReading file " << ligne << "\t";
    this->nbFile++;
    long nbLigneFile = 0;
    short nbOfConcernedIndividual = 0;
    PLOG(plog::debug) << "### size of trnAssociation : " << this->trnAssociation.size(); 
    while (getline(cnvStream, ligneCNV)) {
      PLOG(plog::debug) << "### NEW LINE ###";
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
          PLOG(plog::error) << "\tParsing VCF line : " << ligneCNV << " failed, passing line"; 
          continue; 
        }

        // TRN count
        if ((res[3] == "TRN") || (res[3] == "BND")) {
          PLOG(plog::debug) << "\tEncountering a TRN"; 
          if (this->getFormat() == "VCF"){
            this->parseVCFLineTRN(ligneCNV);
          }
          continue;
        }


        // pass if not del or dup 
        if ((res[3] != "DEL") and (res[3] != "DUP") and (res[3] != "INV")) {
          PLOG(plog::debug) << "\tPassing VCF line : not DEL nor DUP nor INV, passing line : " << res[3]; 
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
      int value = string_to_int(res[4]);
      if (value > 5) {
        value = 5;
      }
      if (res[4] == ".") {
        value = 6;
      }

      PLOG(plog::verbose) << "\tCnv single value for this CNV is " << value;
      
      levelValues[0] = string_to_int(parseOnSep(res[6], ",")[0]);
      levelValues[1] = string_to_int(parseOnSep(res[6], ",")[1]);
      levelValues[2] = string_to_int(parseOnSep(res[6], ",")[2]);
      levelValues[3] = string_to_int(parseOnSep(res[6], ",")[3]);
      levelValues[4] = string_to_int(parseOnSep(res[6], ",")[4]);
      levelValues[5] = string_to_int(parseOnSep(res[6], ",")[5]);
      levelValues[6] = string_to_int(parseOnSep(res[6], ",")[6]); // inversion level

      
      if (value != -1) {
        levelValues[value] = 1;
      }


      for (int cn = 0 ; cn <= 6 ; cn ++) {
        int count = 0; 
        PLOG(plog::debug) << "\tCnv values for this CNV ; cn = " <<  cn  << ", counts = " << levelValues[cn];
        count = levelValues[cn]; 
        if (count == 0) {
          PLOG(plog::debug) << "\t\tnothing to insert";
          continue; 
        }
        
        PLOG(plog::debug) << "\t\twill insert : " << chromosome << ":" << start << "-" << end << " ; cnv : " << cn;
        
        // fill empty map if chr is not existing
        if (!(this->breakpoints.count(chromosome) > 0)) {
          PLOG(plog::debug) << "\t\t\tCreating breakpoint map for this chromosome"; 
          unordered_map<unsigned int, map<long, short> > tempMap;
          this->breakpoints[chromosome] = tempMap;
          map<long, short> tempList;
          this->breakpoints[chromosome][0] = tempList;
          this->breakpoints[chromosome][1] = tempList;
          this->breakpoints[chromosome][2] = tempList;
          this->breakpoints[chromosome][3] = tempList;
          this->breakpoints[chromosome][4] = tempList;
          this->breakpoints[chromosome][5] = tempList;
          this->breakpoints[chromosome][6] = tempList; // inversion level 
        }

        // look for the start / end values
        // if the map is empty do not try to browse it, just insert the start and end values and treat the next line. 
        if (this->breakpoints[chromosome][cn].empty()) {
          PLOG(plog::debug) << "\t\t\tMap was empty : so just inserting start & end";
          this->breakpoints[chromosome][cn][start] = 1;
          this->breakpoints[chromosome][cn][end] = 0;
          continue; 
        }

        // insert start and end values in the sorted map
        short lastCount = 0;
        map<long, short>::iterator it, inserted_it, it_beforestart, it_beforeend, it_afterstart, it_afterend;

        // need to get old values
        PLOG(plog::debug) << "Looking for breakpoint around start and end : ";
        it_beforestart = this->breakpoints[chromosome][cn].lower_bound(start);
        it_beforeend = this->breakpoints[chromosome][cn].lower_bound(end);

        it_afterstart = this->breakpoints[chromosome][cn].upper_bound(start);
        it_afterend = this->breakpoints[chromosome][cn].upper_bound(end);

        if (it_beforestart != this->breakpoints[chromosome][cn].end()) {
          //it_beforestart --;
          PLOG(plog::debug) << "\t\tBreakpoint before start is " << it_beforestart->first << ":" << it_beforestart->second;
        } else {
          PLOG(plog::debug) << "\t\tBreakpoint before start is after the current end of the map";
        }
        if (it_beforeend != this->breakpoints[chromosome][cn].end()) {
          // it_beforeend --;
          PLOG(plog::debug) << "\t\tBreakpoint before end is " << it_beforeend->first << ":" << it_beforeend->second;
        } else {
          PLOG(plog::debug) << "\t\tBreakpoint before end is after the current end of the map";
        }
        if (it_afterstart != this->breakpoints[chromosome][cn].end()) {
          PLOG(plog::debug) << "\t\tBreakpoint after start is " << it_afterstart->first << ":" << it_afterstart->second;
        } else {
          PLOG(plog::debug) << "\t\tBreakpoint after start is after the current end of the map";
        }
        if (it_afterend != this->breakpoints[chromosome][cn].end()) {
          PLOG(plog::debug) << "\t\tBreakpoint after end is " << it_afterend->first << ":" << it_afterend->second;
        } else {
          PLOG(plog::debug) << "\t\tBreakpoint after end is after the current end of the map";
        }

        // manage begin of the map
        if (((it_beforestart == this->breakpoints[chromosome][cn].begin()) && (start != 1)) || (it_beforestart == this->breakpoints[chromosome][cn].end())){
          lastCount = 0;
        } else {
          lastCount =  it_beforestart->second;
        }

        // insert start point 
        this->breakpoints[chromosome][cn].insert_or_assign(start, lastCount + count);
        PLOG(plog::debug) << "\t\tStart inserted " << start << ":" << lastCount + count << " at cn level " << cn;

        // get last value of the interval & manage end of the map
        if (it_beforeend == this->breakpoints[chromosome][cn].end()) {
          lastCount = 0;
        } else {
          lastCount =  it_beforeend->second;
        }

        // modify all value until end
        for (it = it_afterstart ; it != it_afterend ; ++ it) {
          if (it != this->breakpoints[chromosome][cn].end()) {
            this->breakpoints[chromosome][cn][it->first] += count;
            PLOG(plog::debug) << "\t\tChanging breakpoints " << it->first << ":" << this->breakpoints[chromosome][cn][it->first] - count << " to " << this->breakpoints[chromosome][cn][it->first] << " at cn level " << cn;
          }
        }

        // insert the end 
        this->breakpoints[chromosome][value].insert_or_assign(end, lastCount);
        PLOG(plog::debug) << "\t\tEnd inserted " << end << ":" << lastCount << " at cn level " << cn;
        PLOG(plog::debug) << "\t\tSize of breakpoints at chr : " << chromosome << " and value " << cn << " : " << this->breakpoints[chromosome][cn].size();

        // a count for large files to be sure that everything went well
        if ((nbLigneFile % 10000) == 0) {
          PLOG(plog::info) << "\t" << nbLigneFile << " events detected, still in progress";
        }
      }
    }
    PLOG(plog::info) << " with " << nbLigneFile << " events detected ";
  }
  PLOG(plog::info) << "Ended with " << this->getNbFile() << " files";
  PLOG(plog::verbose) << "Leaving cnvCompare::getDataFast ";
}


/**
 * @brief Getter for the number of input files
 * @param none
 * @return integer value : the number of file from the instance of the cnvCompare class
 **/
short int cnvCompare::getNbFile() { 
  PLOG(plog::verbose) << "Entering cnvCompare::getNbFile ";
  return this->nbFile; 
  PLOG(plog::verbose) << "Leaving cnvCompare::getNbFile ";
}


/**
 * @brief Getter for the name of the input file list
 * @param none
 * @return string value : the name of the input file list
 **/
string cnvCompare::getInputFile() {
  PLOG(plog::verbose) << "Entering cnvCompare::getInputFile ";
  return this->inputFile;
  PLOG(plog::verbose) << "Leaving cnvCompare::getInputFile ";
}


/**
 * @brief Getter for the name of the control file list
 * @param none
 * @return string value : the name of the control file list
 **/
string cnvCompare::getControlFile() {
  PLOG(plog::verbose) << "Entering cnvCompare::getControlFile ";
  return this->controlFile;
  PLOG(plog::verbose) << "Leaving cnvCompare::getControlFile ";
}


/**
 * @brief Method used to fill the file map
 * @param incFile : string value : the path to the file containing a list of  files 
 * @param status : string value : precising if the files are controls or inputs
 * @return int value : a return code : 0 if ok, > 0 if not ok
 **/
int cnvCompare::fillMap(string incFile, string status) {
  PLOG(plog::verbose) << "Entering cnvCompare::fillMap ";
  ifstream allStream(incFile.c_str());
  string ligne;
  if (allStream) {
    while (getline(allStream, ligne)) {
      ifstream cnvStream(ligne.c_str());
      PLOG(plog::info) << "\tadding file " << ligne << " as " << status;
      if (!IsFileReadable(ligne)) {
        PLOG(plog::warning) << "File provided as " << status << " : " << ligne << " is not accessible : passsing file";
        continue;
      }
      this->nbFile++;
      this->fileMap[ligne] = status;
    }
  }
  
  PLOG(plog::verbose) << "Leaving cnvCompare::fillMap ";
  return this->nbFile;
}


/**
 * @brief Getter for filter size value
 * @param none
 * @return int value : the filter size value
 **/
int cnvCompare::getFilterSize() {
  PLOG(plog::verbose) << "Entering cnvCompare::getFilterSize ";
  return this->filterSize;
  PLOG(plog::verbose) << "Leaving cnvCompare::getFilterSize ";
}


/**
 * @brief Setter for the file format
 * @param incVCFChoice A boolean indicating if the files are in VCF format
 * @param incBEDChoice A boolean indicating if the files are in BED format
 * @return none
 **/
void cnvCompare::setFormat(bool incVCFChoice, bool incBEDChoice) {
  PLOG(plog::verbose) << "Entering cnvCompare::setFormat ";
  this->useVCFFormat = incVCFChoice;
  this->useBEDFormat = incBEDChoice;
  PLOG(plog::verbose) << "Leaving cnvCompare::setFormat ";
}


/**
 * @brief Getter for the file formate
 * @param none
 * @return string value : the file format
 **/
string cnvCompare::getFormat() {
  PLOG(plog::verbose) << "Entering cnvCompare::getFormat ";
  if (this->useVCFFormat) {
    PLOG(plog::verbose) << "Leaving cnvCompare::getFormat ";
    return "VCF";
  } else {
    PLOG(plog::verbose) << "Leaving cnvCompare::getFormat ";
    return "BED";
  }
  PLOG(plog::verbose) << "Leaving cnvCompare::getFormat ";
}


/**
 * @brief Method used to populate the chromosome map
 * @param none
 * @return int value : a return code : 0 if ok, > 0 if not ok
 **/
int cnvCompare::populateChr() {
  PLOG(plog::verbose) << "Entering cnvCompare::populateChr ";
  this->chromosomeMap.insert(this->chromosomeMap.end(), {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"});
  PLOG(plog::verbose) << "Leaving cnvCompare::populateChr ";
  return 0; 
}


/**
 * @brief Setter used to define the suffix to use
 * @param incSuffix : string defining the suffix value
 * @return none
 **/
void cnvCompare::setSuffix(string incSuffix) {
  PLOG(plog::verbose) << "Entering cnvCompare::setSuffix ";
  this->suffix = incSuffix;
  PLOG(plog::verbose) << "Leaving cnvCompare::setSuffix ";
}


/**
 * @brief getter used to get the suffix to use
 * @param none
 * @return string value : defining the suffix value
 **/
string cnvCompare::getSuffix() {
  PLOG(plog::verbose) << "Entering cnvCompare::getSuffix ";
  PLOG(plog::verbose) << "Leaving cnvCompare::getSuffix ";
  return this->suffix; 
}


/**
 * @brief Setter used to define the dict file to use to populate the chromosome map
 * @param incDictFile : string defining the dict file path
 * @return none
 **/
void cnvCompare::setDictFile(string incDictFile) {
  PLOG(plog::verbose) << "Entering cnvCompare::setDictFile ";
  this->dictFile = incDictFile;
  PLOG(plog::verbose) << "Leaving cnvCompare::setDictFile ";
}


/**
 * @brief getter used to get the dict file to use
 * @param none
 * @return string value defining the dict file
 **/
string cnvCompare::getDictFile() {
  PLOG(plog::verbose) << "Entering cnvCompare::getDictFile ";
  PLOG(plog::verbose) << "Leaving cnvCompare::getDictFile ";
  return this->dictFile; 
}


/**
 * @brief Method used to parse a dictionnary file and to fill the apporpriate container. 
 * @param incDictFile : string defining path to a dictionnary file
 * @return none
 **/
void cnvCompare::parseDictFile(string incDictFile) {
  PLOG(plog::verbose) << "Entering cnvCompare::parseDictFile ";
  string ligne;
  string mot;
  string header = "@SQ";

  // Parsing dict file
  PLOG(plog::info) << "Reading Dict file " << incDictFile;
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
      PLOG(plog::info) << "adding " << chrName << " to the chromosome map";
      this->chromosomeMap.push_back(chrName);
    }
  }
  PLOG(plog::verbose) << "Leaving cnvCompare::parseDictFile ";
}

/**
 * @brief Setter used to define if a dict file is present
 * @param incBool : boolean indicating if a dict file is present
 * @return none
 **/
void cnvCompare::setHasDict(bool incBool) {
  PLOG(plog::verbose) << "Entering cnvCompare::setHasDict ";
  this->hasDict = incBool;
  PLOG(plog::verbose) << "Leaving cnvCompare::setHasDict ";
}


/**
 * @brief Getter for number of individuals (different from number of files)
 * @param none
 * @return short value : the number of individual in memory
 **/
short int cnvCompare::getNbIndividual() {
  PLOG(plog::verbose) << "Entering cnvCompare::getNbIndividual ";
  PLOG(plog::verbose) << "Leaving cnvCompare::getNbIndividual ";
  return this->nbIndividual;
}

/**
 * @brief Method used to get the number of individuals on a header line for a VCF
 * @param incLine : A string containing the header line
 * @return none
 **/
void cnvCompare::watchHeader(string incLine) {
  PLOG(plog::verbose) << "Entering cnvCompare::watchHeader ";
  string goodHeader = "#CHR"; 
  short initIndiv = this->getNbIndividual();
  if (incLine.find(goodHeader) == 0) {
    vector<string> lineTable = parseOnSep(incLine, "\t");
    for (long unsigned int n = 9 ; n < lineTable.size() ; n ++) {
      this->nbIndividual += 1;
    }
    PLOG(plog::info) << "Added " << this->getNbIndividual() - initIndiv << " individuals to the list, total is now : " << this->getNbIndividual();
  }

  PLOG(plog::verbose) << "Leaving cnvCompare::watchHeader";
}


/**
 * @brief Method used to insert a breakpoint into the dedicated container
 * @brief Will modify the breakpoints if existing
 * @param incMap : A map containing the breakpoint informations
 * @return integer value : the index of the breakpoint for the future association
 **/
 int cnvCompare::insertTrnBreakpoints(string chromosome, unsigned int position) {
  PLOG(plog::verbose) << "Entering cnvCompare::insertTrnBreakpoints ";
  
  // compute next index if needed 
  unordered_map<string, map<unsigned int, int> >::iterator myIterA;
  map<unsigned int, int>::iterator myIterB;
  map<unsigned int, int> tmpMap;
  int index = this->getNextIndex();

  // fill empty map if chr is not existing
  if (!(this->trnbreakpoints.count(chromosome) > 0)) {
    PLOG(plog::info) << "\tCreating TRN breakpoint map for this chromosome " << chromosome; 
    map<unsigned int, int> tempMap;
    this->trnbreakpoints[chromosome] = tempMap;
  }

  // look for the start / end values
  // if the map is empty do not try to browse it, just insert the start and end values and treat the next line. 
  if (this->trnbreakpoints[chromosome].empty()) {
    PLOG(plog::verbose) << "\t\t\tMap was empty : so just inserting Position ";
    this->trnbreakpoints[chromosome][position] = index;
    this->setNextIndex(index + 1);
    PLOG(plog::verbose) << "Leaving cnvCompare::insertTrnBreakpoints";
    return index; 
  }

  // insert start and end values in the sorted map
  map<unsigned int, int>::iterator it_beforeposition, it_afterposition;

  // need to get values
  it_beforeposition = this->trnbreakpoints[chromosome].lower_bound(position);
  it_afterposition = this->trnbreakpoints[chromosome].upper_bound(position);


  if (it_beforeposition != this->trnbreakpoints[chromosome].begin()) {
    PLOG(plog::verbose) << "\tBreakpoint before position is " << it_beforeposition->first << " index :" << it_beforeposition->second;
    if ((position >= (it_beforeposition->first - 200)) && (position <= (it_beforeposition->first + 200))) {
      PLOG(plog::verbose) << "\t\tBreakpoint is the same";
      PLOG(plog::verbose) << "Leaving cnvCompare::insertTrnBreakpoints";
      return it_beforeposition->second;
    }
  } else {
    PLOG(plog::verbose) << "\tBreakpoint before position is before the current begin of the map";
  }

  if (it_afterposition != this->trnbreakpoints[chromosome].end()) {
    PLOG(plog::verbose) << "\tBreakpoint after position is " << it_afterposition->first << " index :" << it_afterposition->second;
    if ((position <= (it_afterposition->first + 200)) && (position >= (it_afterposition->first - 200))) {
      PLOG(plog::verbose) << "\t\tBreakpoint is the same";
      PLOG(plog::verbose) << "Leaving cnvCompare::insertTrnBreakpoints";
      return it_afterposition->second;
    }
  } else {
    PLOG(plog::verbose) << "\tBreakpoint after position is after the current end of the map";
  }


  PLOG(plog::verbose) << "\tBreakpoint wasn't found, inserting it";
  this->trnbreakpoints[chromosome][position] = index;
  this->setNextIndex(index + 1);
  PLOG(plog::verbose) << "Leaving cnvCompare::insertTrnBreakpoints";
  return index; 


  PLOG(plog::verbose) << "Leaving cnvCompare::insertTrnBreakpoints";
  return 0;
 }

/**
 * @brief Method used to get an index from a position in the trn breakpoint list
 * @brief does not insert or modify the breakpoints
 * @param chromosome : a string containeing chromosome
 * @param position : unsigned int containing the position on chromosome
 * @return integer value : the index of the breakpoint
 **/
int cnvCompare::getTRNIndex(string chromosome, unsigned int position) {
  PLOG(plog::verbose) << "Entering cnvCompare::getTRNIndex ";

  // itor to get the values 
  map<unsigned int, int>::iterator it_beforeposition, it_afterposition;

  // need to get values
  it_beforeposition = this->trnbreakpoints[chromosome].lower_bound(position);
  it_afterposition = this->trnbreakpoints[chromosome].upper_bound(position);

  if (it_beforeposition != this->trnbreakpoints[chromosome].begin()) {
    if ((position >= (it_beforeposition->first - 200)) && (position <= (it_beforeposition->first + 200))) {
      PLOG(plog::verbose) << "Leaving cnvCompare::insertTrnBreakpoints";
      PLOG(plog::debug) << "\t\tFound (before) : " << chromosome << ":" << it_beforeposition->first << " => " << it_beforeposition->second;
      return it_beforeposition->second;
    }
  } 

  if (it_afterposition != this->trnbreakpoints[chromosome].end()) {
    if ((position <= (it_afterposition->first + 200)) && (position >= (it_afterposition->first - 200))) {
      PLOG(plog::verbose) << "Leaving cnvCompare::insertTrnBreakpoints";
      PLOG(plog::debug) << "\t\tFound (after) : " << chromosome << ":" << it_afterposition->first << " => " << it_afterposition->second;
      return it_afterposition->second;
    }
  }
  PLOG(plog::debug) << "\t\tBreakpoint not found in the map => 1";
  PLOG(plog::verbose) << "Leaving cnvCompare::getTRNIndex";
  return 1;
}

/**
 * @brief Method used to get the next index available
 * @param None
 * @return integer value : the next index available
 **/
int cnvCompare::getNextIndex() {
  PLOG(plog::verbose) << "Entering cnvCompare::getTRNIndex";
  return this->nextIndex;
  PLOG(plog::verbose) << "Leaving cnvCompare::getNextIndex";
}

/**
 * @brief Method used to set the next index available
 * @param Integer being the next index 
 * @return void
 **/
 void cnvCompare::setNextIndex(int incValue) {
  PLOG(plog::verbose) << "Entering cnvCompare::setTRNIndex";
  this->nextIndex = incValue;
  PLOG(plog::verbose) << "Leaving cnvCompare::setNextIndex";
}

/**
 * @brief Method used to infer copy number with a GT field
 * @param String the GT Field
 * @return Int the value of the copynumber 
 **/
int cnvCompare::inferCNfromGT(string GT, string svType) {
  PLOG(plog::verbose) << "Entering cnvCompare::inferCNfromGT";

  if ((GT == "./.") || (GT == ".")){
    PLOG(plog::debug) << "GT field is ./. : CN = 2";
    return 2;
  }
  
  int al1 = string_to_int(parseOnSep(GT, "/")[0]);
  int al2 = string_to_int(parseOnSep(GT, "/")[1]);
  PLOG(plog::debug ) << "Al 1 = " << al1 << " ; Al 2 = " << al2;
  if ((al1 == 0) && (al2 == 0)) {
    return 2;
  }

  if (svType == "DEL") {
    if (al1 != al2) {
      return 1;
    } else {
      return 0;
    }
  }
  if (svType == "DUP") {
    if (al1 != al2) {
      return 3;
    } else {
      return 4;
    }
  }

  
  PLOG(plog::verbose) << "Leaving cnvCompare::inferCNfromGT";
  return 2;
}