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
#include <iomanip>
#include <unordered_map>
//#include <ranges>

// Boost 
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/expressions.hpp>


// utils
#include "utils.h"

// Classes
#include "cnvCompare.h"

using namespace std;
using namespace boost;
namespace logging = boost::log;

/**
 * @brief default constructor (useless)
 * @param none
 * @return none
 **/
cnvCompare::cnvCompare()
{
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::cnvCompare default constructor" << endl;
  //
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
  BOOST_LOG_TRIVIAL(trace) <<  "sizeof( char ) = " << sizeof(char);
  BOOST_LOG_TRIVIAL(trace) <<  "sizeof( string ) = " << sizeof(string);
  BOOST_LOG_TRIVIAL(trace) <<  "sizeof( int32_t ) = " << sizeof(int32_t);
  BOOST_LOG_TRIVIAL(trace) <<  "sizeof( int ) = " << sizeof(int);

  this->inputFile = iF;
  this->nbThread = nT;
  this->useControls = false;
  this->filterSize = s;
  this->suffix = "count";
  int n = this->fillMap(iF, "sample");
  BOOST_LOG_TRIVIAL(info) << n << " files added to sample list" << endl;
  int ret = this->populateChr();
  if (ret != 0) {
    BOOST_LOG_TRIVIAL(error) << "Chromosome population couldn't be filled in" << endl;
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
  this->inputFile = iF;
  this->controlFile = cF;
  this->nbThread = nT;
  this->useControls = true;
  this->filterSize = s;
  this->suffix = "count";
  int n = this->fillMap(iF, "sample");
  BOOST_LOG_TRIVIAL(info) << n << " files added to sample list" << endl;
  int m = this->fillMap(cF, "control");
  BOOST_LOG_TRIVIAL(info) << m << " files added to control list" << endl;
  int ret = this->populateChr();
  if (ret != 0) {
    BOOST_LOG_TRIVIAL(error) << "Chromosome population couldn't be filled in" << endl;
  }
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::cnvCompare (string, string, int, int) Ctor" << endl;
}


/**
 * @brief Main loop used to get data, and compute counts on the whole genome ; need huge amounts of RAM
 * @param none
 * @return none
 **/
void cnvCompare::mainLoop()
{
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::mainLoop " << endl;
  this->getData();
  this->computeCounts();
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::mainLoop " << endl;
}

/**
 * @brief Alt loop used to get data, and compute counts chr by chr
 * @param none
 * @return none
 **/
void cnvCompare::altLoop()
{
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::altLoop " << endl;
  for (auto &a : this->chromosomeMap)
  {
    this->getDataWhole(a);
    this->computeChrCounts(a);
    this->cleanData();
  }
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::altLoop " << endl;
}

/**
 * @brief Method used to emptying data for the chromosomes
 * @param none
 * @return none
 **/
void cnvCompare::cleanData() { 
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::cleanData " << endl;
  this->dataByChr.clear();
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::cleanData " << endl; 
}

/**
 * @brief Method used to collect data from input files
 * @param incChr the chr to filter
 * @return none
 * @todo Need to clarify which method is used chr by chr, actually the 2 loops are using the chr by chr algo... 
 **/
void cnvCompare::getDataWhole(string incChr)
{
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getDataWhole " << endl;
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
  for (myIterA = this->fileMap.begin(); myIterA != this->fileMap.end(); myIterA++)
  {
    ligne = myIterA->first;
    ifstream cnvStream(ligne.c_str());
    BOOST_LOG_TRIVIAL(info) << "\tReading file " << ligne << "\t" << endl;
    this->nbFile++;
    long nbLigneFile = 0;
    while (getline(cnvStream, ligneCNV))
    {

      // need to deal with header "#"
      if (ligneCNV.find(header) == 0)
      {
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
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getDataWhole " << endl;
}


/**
 * @brief Method used to commpute counts from in memory data chr by chr
 * @param incChr the chr to filter
 * @return none
 * @todo Need to clarify which method is used chr by chr, actually the 2 loops are using the chr by chr algo... 
 **/
void cnvCompare::computeChrCounts(string incChr)
{
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::computeChrCounts " << endl;
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

  // tsv parsing
  map<string, string>::iterator myIterA;
  for (myIterA = this->fileMap.begin(); myIterA != this->fileMap.end(); myIterA++)
  {
    // struct timeval tbegin, tend;
    ligne = myIterA->first;
    string s = myIterA->second;
    if (s == "control")
    {
      continue;
    }
    ifstream cnvStream(ligne.c_str());
    BOOST_LOG_TRIVIAL(info) << "\tReading file " << ligne << endl;
    // if "merged" isn't present in filename, renaming failed and rewrite the original file !
    // TO CORRECT
    string outFileName = pyReplace(ligne, "merged", "count");
    BOOST_LOG_TRIVIAL(info) << "\t\tWriting in file " << outFileName << endl;
    ofstream outStream;

    if (incChr == "chr1")
    {
      outStream.open(outFileName.c_str(), ios::out);
    }
    else
    {
      outStream.open(outFileName.c_str(), ios::out | ios::app);
    }


    // cerr << "The format is : " << this->getFormat() << endl; 

    while (getline(cnvStream, ligneCNV))
    {
      // need to deal with header "#"
      if (ligneCNV.find(header) == 0)
      {
        if ((this->getFormat() == "VCF") && (incChr == "chr1"))
        {
          outStream << ligneCNV << endl;
        }
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
      }

      // type conversion
      chromosome = res[0];
      s_type = res[3];
      long start = string_to_int(res[1]);
      long end = string_to_int(res[2]);
      unsigned int value = string_to_int(res[4]);
      if (value > 5)
      {
        value = 5;
      }
      // skip if not the good chr
      if (chromosome != incChr)
      {
        continue;
      }

      // pass if not del or dup
      BOOST_LOG_TRIVIAL(debug) << "s_type = " << s_type << endl;
      if ((s_type != "DUP") && (s_type != "DEL"))
      {
        outStream << ligneCNV << endl;
        continue;
      }

      // counts 
      double total = 0;
      for (long i = start; i <= end; i++)
      {
        total += (double)(this->dataByChr[value][i]);
      }
      double mean = total / (double)((end-start)+1);

      // need to adapt the output according to the choosen format
      if (this->getFormat() == "BED")
      {
        outStream << chromosome << "\t" << start << "\t" << end << "\t";
        if (value > 2)
        {
          outStream << "DUP\t";
        }
        else
        {
          outStream << "DEL\t";
        }
        outStream << value << "\t" << mean << "/" << this->getNbFile() << endl;
      }
      else
      {
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

        while (getline(issLigne, mot, '\t'))
        {
          switch (i)
          {
          case 0:
            outStream << mot;
            break;
          case 7:
            outStream << "\t";
            info = mot;
            issInfo.str(info);
            while (getline(issInfo, infomot, ';'))
            {
              if (infomot.find("SVTYPE=") == 0)
              {
                svtype = parseOnSep(infomot, "=")[1];
                continue;
              }
              if (infomot.find("END=") == 0)
              {
                ciend = parseOnSep(infomot, "=")[1];
                continue;
              }
              if (infomot.find("VALUE=") == 0)
              {
                value = parseOnSep(infomot, "=")[1];
                continue;
              }
              outStream << infomot << ";";
            }
            outStream << "END=" << ciend << ";VALUE=" << value << ";SVTYPE=";
            if (string_to_int(value) > 2)
            {
              outStream << "DUP;";
            }
            else
            {
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
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::computeChrCounts " << endl;
}


/**
 * @brief Method used to parse a line from a BED file
 * @param incLine the line to parse
 * @return vector<string> containing values : chr, start, end, type, value
 * @todo need to check that all the values are present
 **/
// output a vector containing : chr start end type value from a BED line
vector<string> cnvCompare::parseBEDLine(string incLine)
{
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::parseBEDLine " << endl;
  vector<string> output;
  string mot;
  short int i = 0;

  // get information from the line
  istringstream issLigne(incLine);
  while (getline(issLigne, mot, '\t'))
  {
    switch (i)
    {
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
vector<string> cnvCompare::parseVCFLine(string incLine)
{
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::parseVCFLine " << endl;
  vector<string> output;
  map<string, string> temp;
  string mot;
  string infomot;
  string info;
  short int i = 0;
  istringstream issInfo;
  temp["SVTYPE"] = "NONE";

  // get information from the line
  istringstream issLigne(incLine);
  i = 0;
  while (getline(issLigne, mot, '\t'))
  {
    switch (i)
    {
    case 0:
      output.push_back(mot);
      break;
    case 1:
      output.push_back(mot);
      break;
    case 7:
      info = mot;
      issInfo.str(info);
      while (getline(issInfo, infomot, ';'))
      {
        if (infomot.find("SVTYPE=") == 0)
        {
          temp["SVTYPE"] = parseOnSep(infomot, "=")[1];
        }
        if (infomot.find("END=") == 0)
        {
          temp["END"] = parseOnSep(infomot, "=")[1];
        }
        if (infomot.find("VALUE=") == 0)
        {
          temp["VALUE"] = parseOnSep(infomot, "=")[1];
        }
      }
      break;
    default:
      break;
    }
    i++;
  }

  // add data to the vector from the temp map
  output.push_back(temp["END"]);
  output.push_back(temp["SVTYPE"]);
  output.push_back(temp["VALUE"]);
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::parseVCFLine " << endl;
  return output;
}

/**
 * @brief Method used to collect data from input files : chr by chr mode
 * @param none
 * @return none
 * @todo Need to clarify which method is used chr by chr, actually the 2 loops are using the chr by chr algo... 
 **/
void cnvCompare::getData()
{
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getData " << endl;
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
        continue;
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
      if ((end - start) < this->getFilterSize())
      {
        continue;
      }

      if (value > 5)
      {
        value = 5;
      }

      if (!(this->data.count(chromosome) > 0))
      {
        map<unsigned int, map<long, vector<long>>> tempMap;
        this->data[chromosome] = tempMap;
        map<long, vector<long>> tempMap2;
        this->data[chromosome][0] = tempMap2;
        this->data[chromosome][1] = tempMap2;
        this->data[chromosome][2] = tempMap2;
        this->data[chromosome][3] = tempMap2;
        this->data[chromosome][4] = tempMap2;
        this->data[chromosome][5] = tempMap2;
      }
      this->data[chromosome][value][start].push_back(end);
    }
    BOOST_LOG_TRIVIAL(info) << " with " << nbLigneFile << " events detected " << endl;
  }
  BOOST_LOG_TRIVIAL(info) << "Ended with " << this->getNbFile() << " files" << endl;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getData " << endl;
}

/**
 * @brief Method used to commpute counts from in memory data on the whole genome
 * @param none
 * @return none
 * @todo Need to clarify which method is used chr by chr, actually the 2 loops are using the chr by chr algo... 
 **/
void cnvCompare::computeCounts()
{
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::computeCounts " << endl;
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

  map<string, string>::iterator myIterA;
  for (myIterA = this->fileMap.begin(); myIterA != this->fileMap.end(); myIterA++)
  {
    ligne = myIterA->first;
    string s = myIterA->second;
    if (s == "control")
    {
      continue;
    }

    ifstream cnvStream(ligne.c_str());
    BOOST_LOG_TRIVIAL(info) << "\tReading file " << ligne << endl;
    string outFileName = pyReplace(ligne, "merged", "count");
    BOOST_LOG_TRIVIAL(info) << "\t\tWriting in file " << outFileName << endl;
    ofstream outStream;
    outStream.open(outFileName.c_str(), ios::out);
    while (getline(cnvStream, ligneCNV))
    {
      // need to deal with header "#"
      if (ligneCNV.find(header) == 0)
      {
        if (this->getFormat() == "VCF")
        {
          outStream << ligneCNV << endl;
        }
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
      }

      // type conversion
      chromosome = res[0];
      s_type = res[3];

      // pass if not del or dup
      // verifier la condition : a priori les BND et INV passent à travers
      // cerr << "DEBUG : s_type = " << s_type << endl;

      if ((s_type != "DUP") && (s_type != "DEL"))
      {
        outStream << ligneCNV << endl;
        continue;
      }

      long start = string_to_int(res[1]);
      long end = string_to_int(res[2]);
      unsigned int value = string_to_int(res[4]);

      if (value > 5)
      {
        value = 5;
      }
      // counts
      map<long, short int> counts;
      map<long, vector<long>> tmpMap;
      tmpMap = this->data[chromosome][value];
      // init data counts
      for (long j = start; j <= end; j++)
      {
        counts[j] = 0;
      }
      // get pair of iterators on range of intervals including the start point
      auto lb = tmpMap.lower_bound(start);
      auto ub = tmpMap.upper_bound(end);
      // cerr << "\t\t\tGet range of starts from : \n";
      map<long, vector<long>>::iterator myIterC;
      for (myIterC = lb; myIterC != ub; myIterC++)
      {
        vector<long>::iterator myIterD;
        for (myIterD = (myIterC->second).begin();
             myIterD != (myIterC->second).end(); myIterD++)
        {
          // cerr <<  "\t\t\t\t" << myIterC->first << " to " << *myIterD <<
          // endl;
          long i = myIterC->first;
          while (i <= *myIterD)
          {
            if (i > end)
            {
              break;
            }
            counts[i]++;
            i++;
          }
        }
      }
      vector<double> m;
      for (auto &s : counts)
      {
        m.push_back((double)(s.second));
      }
      double mean = moyenne_calculator(m);

      // need to adapt the output according to the choosen format
      if (this->getFormat() == "BED")
      {
        outStream << chromosome << "\t" << start << "\t" << end << "\t";
        if (value > 2)
        {
          outStream << "DUP\t";
        }
        else
        {
          outStream << "DEL\t";
        }
        outStream << value << "\t" << mean << "/" << this->getNbFile() << endl;
      }
      else
      {
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

        while (getline(issLigne, mot, '\t'))
        {
          switch (i)
          {
          case 0:
            outStream << mot;
            break;
          case 7:
            outStream << "\t";
            info = mot;
            issInfo.str(info);
            while (getline(issInfo, infomot, ';'))
            {
              if (infomot.find("SVTYPE=") == 0)
              {
                svtype = parseOnSep(infomot, "=")[1];
                continue;
              }
              if (infomot.find("END=") == 0)
              {
                ciend = parseOnSep(infomot, "=")[1];
                continue;
              }
              if (infomot.find("VALUE=") == 0)
              {
                value = parseOnSep(infomot, "=")[1];
                continue;
              }
              outStream << infomot << ";";
            }
            outStream << "END=" << ciend << ";VALUE=" << value << ";SVTYPE=";
            if (string_to_int(value) > 2)
            {
              outStream << "DUP;";
            }
            else
            {
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
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::computeCounts " << endl;
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
string cnvCompare::getInputFile()
{
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getInputFile " << endl;
  return this->inputFile;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getInputFile " << endl;
}

/**
 * @brief Getter for the name of the control file list
 * @param none
 * @return string value : the name of the control file list
 **/
string cnvCompare::getControlFile()
{
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
int cnvCompare::fillMap(string incFile, string status)
{
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::fillMap " << endl;
  ifstream allStream(incFile.c_str());
  string ligne;
  if (allStream)
  {
    while (getline(allStream, ligne))
    {
      ifstream cnvStream(ligne.c_str());
      BOOST_LOG_TRIVIAL(info) << "\tadding file " << ligne << " as " << status << endl;
      if (!IsFileReadable(ligne))
      {
        BOOST_LOG_TRIVIAL(warning) << "File provided as " << status << " : " << ligne << " is not accessible : passsing file" << endl;
        continue;
      }

      this->nbFile++;
      this->fileMap[ligne] = status;
    }
  }
  return this->nbFile;
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::fillMap " << endl;
}

/**
 * @brief Getter for filter size value
 * @param none
 * @return int value : the filter size value
 **/
int cnvCompare::getFilterSize()
{
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
void cnvCompare::setFormat(bool incVCFChoice, bool incBEDChoice)
{
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
string cnvCompare::getFormat()
{
  BOOST_LOG_TRIVIAL(trace) << "Entering cnvCompare::getFormat " << endl;
  if (this->useVCFFormat)
  {
    return "VCF";
  }
  else
  {
    return "BED";
  }
  BOOST_LOG_TRIVIAL(trace) << "Leaving cnvCompare::getFormat " << endl;
}

/**
 * @brief Method used to populate the chromosome map
 * @param none
 * @return int value : a return code : 0 if ok, > 0 if not ok
 * @todo Need to deal with non human genomes
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