// C++ std libs
#include <iostream>
// #include <iomanip>
#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/time.h>
#include <sys/types.h>
#include <vector>

// utils
#include "utils.h"

// Classes
#include "cnvCompare.h"

using namespace std;

cnvCompare::cnvCompare() {
  //
}

cnvCompare::cnvCompare(string iF, int nT, int s) {
  this->inputFile = iF;
  this->nbThread = nT;
  this->useControls = false;
  this->filterSize = s;
  int n = this->fillMap(iF, "sample");
  cerr << n << " files added to sample list" << endl; 
  this->chromosomeMap.push_back("chr1");
  this->chromosomeMap.push_back("chr2");
  this->chromosomeMap.push_back("chr3");
  this->chromosomeMap.push_back("chr4");
  this->chromosomeMap.push_back("chr5");
  this->chromosomeMap.push_back("chr6");
  this->chromosomeMap.push_back("chr7");
  this->chromosomeMap.push_back("chr8");
  this->chromosomeMap.push_back("chr9");
  this->chromosomeMap.push_back("chr10");
  this->chromosomeMap.push_back("chr11");
  this->chromosomeMap.push_back("chr12");
  this->chromosomeMap.push_back("chr13");
  this->chromosomeMap.push_back("chr14");
  this->chromosomeMap.push_back("chr15");
  this->chromosomeMap.push_back("chr16");
  this->chromosomeMap.push_back("chr17");
  this->chromosomeMap.push_back("chr18");
  this->chromosomeMap.push_back("chr19");
  this->chromosomeMap.push_back("chr20");
  this->chromosomeMap.push_back("chr21");
  this->chromosomeMap.push_back("chr22");
  this->chromosomeMap.push_back("chrX");
  this->chromosomeMap.push_back("chrY");
}

cnvCompare::cnvCompare(string iF, string cF, int nT, int s) {
  this->inputFile = iF;
  this->controlFile = cF;
  this->nbThread = nT;
  this->useControls = true;
  this->filterSize = s;
  int n = this->fillMap(iF, "sample");
  cerr << n << " files added to sample list" << endl; 
  int m = this->fillMap(cF, "control");
  cerr << m << " files added to control list" << endl; 

  this->chromosomeMap.push_back("chr1");
  this->chromosomeMap.push_back("chr2");
  this->chromosomeMap.push_back("chr3");
  this->chromosomeMap.push_back("chr4");
  this->chromosomeMap.push_back("chr5");
  this->chromosomeMap.push_back("chr6");
  this->chromosomeMap.push_back("chr7");
  this->chromosomeMap.push_back("chr8");
  this->chromosomeMap.push_back("chr9");
  this->chromosomeMap.push_back("chr10");
  this->chromosomeMap.push_back("chr11");
  this->chromosomeMap.push_back("chr12");
  this->chromosomeMap.push_back("chr13");
  this->chromosomeMap.push_back("chr14");
  this->chromosomeMap.push_back("chr15");
  this->chromosomeMap.push_back("chr16");
  this->chromosomeMap.push_back("chr17");
  this->chromosomeMap.push_back("chr18");
  this->chromosomeMap.push_back("chr19");
  this->chromosomeMap.push_back("chr20");
  this->chromosomeMap.push_back("chr21");
  this->chromosomeMap.push_back("chr22");
  this->chromosomeMap.push_back("chrX");
  this->chromosomeMap.push_back("chrY");
}

void cnvCompare::mainLoop() {
  this->getData();
  this->computeCounts();
}




// Alt loop : store all data for 1 chromosome and perform counts on it
// need a huge amount of RAM.
void cnvCompare::altLoop() {
  for (auto &a : this->chromosomeMap) {
    this->getDataWhole(a);
    this->computeChrCounts(a);
    this->cleanData();
  }
}

void cnvCompare::cleanData() { this->dataByChr.clear(); }

void cnvCompare::getDataWhole(string incChr) {
  cerr << "Gathering data for chr " << incChr << endl;
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

  // prefill the container
  map<long, short> tempMap2;
  this->dataByChr[0] = tempMap2;
  this->dataByChr[1] = tempMap2;
  this->dataByChr[2] = tempMap2;
  this->dataByChr[3] = tempMap2;
  this->dataByChr[4] = tempMap2;
  this->dataByChr[5] = tempMap2;

  this->nbFile = 0;

  // tsv parsing
  map <string, string>::iterator myIterA;
  for (myIterA = this->fileMap.begin(); myIterA != this->fileMap.end(); myIterA++) {
    ligne = myIterA->first;
    ifstream cnvStream(ligne.c_str());
    cerr << "\tReading file " << ligne << "\t";
    this->nbFile++;
    long nbLigneFile = 0;
    while (getline(cnvStream, ligneCNV)) {

      // need to deal with header "#"
      if (ligneCNV.find(header) == 0) {
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
      
      // Pass the line if not the asked chromosome
      if (chromosome != incChr) {
        continue;
      }
      nbLigneFile++;
      s_type = res[3];
      long start = string_to_int(res[1]);
      long end = string_to_int(res[2]);
      unsigned int value = string_to_int(res[4]);

      // size filter 
      if ((end - start) < this->getFilterSize()) {
        continue;
      }

      // thresholding the dups
      if (value > 5) {
        value = 5;
      }

      for (long i = start; i <= end; i++) {
        if (not((this->dataByChr[value]).count(i) > 0)) {
          this->dataByChr[value][i] = 0;
        }
        this->dataByChr[value][i]++;
      }
    }
    cerr << " with " << nbLigneFile << " events detected " << endl;
  }
  cerr << "Ended with " << this->getNbFile() << " files" << endl;
}

void cnvCompare::computeChrCounts(string incChr) {
  std::cout.precision(3);
  cerr << "Computing counts for chr " << incChr << endl;
  struct timeval tbegin, tend;
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
  map <string, string>::iterator myIterA;
  for (myIterA = this->fileMap.begin(); myIterA != this->fileMap.end(); myIterA++) {
    ligne = myIterA->first;
    string s = myIterA->second; 
    if (s == "control") {
      continue;
    }
    ifstream cnvStream(ligne.c_str());
    cerr << "\tReading file " << ligne << endl;
    // if "merged" isn't present in filename, renaming failed and rewrite the original file !
    // TO CORRECT 
    string outFileName = pyReplace(ligne, "merged", "count");
    cerr << "\t\tWriting in file " << outFileName << endl;
    ofstream outStream;

    if (incChr == "chr1") {
      outStream.open(outFileName.c_str(), ios::out);
    } else {
      outStream.open(outFileName.c_str(), ios::out | ios::app);
    }

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


      // counts
      vector<double> m;
      for (long i = start; i <= end; i++) {
        if ((this->dataByChr[value]).count(i) > 0) {
          m.push_back((double)(this->dataByChr[value][i]));
        } else {
          m.push_back(0.0);
        }
      }

      double mean = moyenne_calculator(m);
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
            while (getline(issInfo, infomot , ';')) {
              if (infomot.find("SVTYPE=") == 0) {
                svtype = parseOnSep(infomot, "=")[1];
                continue;
              }
              if (infomot.find("CIEND=") == 0) {
                ciend = parseOnSep(infomot, "=")[1];
                continue;
              }
              if (infomot.find("VALUE=") == 0) {
                value = parseOnSep(infomot, "=")[1];
                continue;
              }
              outStream << infomot << ";";
            }
            outStream << "CIEND=" << ciend << ";VALUE=" << value << ";SVTYPE="; 
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
}


// output a vector containing : chr start end type value from a BED line
vector<string> cnvCompare::parseBEDLine(string incLine) {
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
}


// output a vector containing : chr start end type value from a VCF line
vector<string> cnvCompare::parseVCFLine(string incLine) {
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
  while (getline(issLigne, mot, '\t')) {
    switch (i) {
    case 0:
      output.push_back(mot);
      break;
    case 1:
      output.push_back(mot);
      break;
    case 7:
      info = mot; 
      issInfo.str(info);
      while (getline(issInfo, infomot , ';')) {
        if (infomot.find("SVTYPE=") == 0) {
          temp["SVTYPE"] = parseOnSep(infomot, "=")[1];
        }
        if (infomot.find("CIEND=") == 0) {
          temp["CIEND"] = parseOnSep(infomot, "=")[1];
        }
        if (infomot.find("VALUE=") == 0) {
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
  output.push_back(temp["CIEND"]);
  output.push_back(temp["SVTYPE"]);
  output.push_back(temp["VALUE"]);

  return output; 
}





void cnvCompare::getData() {
  cerr << "Gathering data" << endl;
  struct timeval tbegin, tend;
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
  map <string, string>::iterator myIterA;
  for (myIterA = this->fileMap.begin(); myIterA != this->fileMap.end(); myIterA++) {
    ligne = myIterA->first;
    ifstream cnvStream(ligne.c_str());
    cerr << "\tReading file " << ligne << "\t";
    this->nbFile++;
    long nbLigneFile = 0;
    while (getline(cnvStream, ligneCNV)) {
      // need to deal with header "#"
      if (ligneCNV.find(header) == 0) {
        continue;
      }
      nbLigneFile++;
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

      // size filter 
      if ((end - start) < this->getFilterSize()) {
        continue;
      }


      if (value > 5) {
        value = 5;
      }

      if (!(this->data.count(chromosome) > 0)) {
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
    cerr << " with " << nbLigneFile << " events detected " << endl;
  }
  cerr << "Ended with " << this->getNbFile() << " files" << endl;
}

void cnvCompare::computeCounts() {
  cerr << "Computing counts" << endl;
  struct timeval tbegin, tend;
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

  map <string, string>::iterator myIterA;
  for (myIterA = this->fileMap.begin(); myIterA != this->fileMap.end(); myIterA++) {
    ligne = myIterA->first;
    string s = myIterA->second; 
    if (s == "control") {
      continue;
    }

    ifstream cnvStream(ligne.c_str());
    cerr << "\tReading file " << ligne << endl;
    string outFileName = pyReplace(ligne, "merged", "count");
    cerr << "\t\tWriting in file " << outFileName << endl;
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
      long start = string_to_int(res[1]);
      long end = string_to_int(res[2]);
      unsigned int value = string_to_int(res[4]);
      
      if (value > 5) {
        value = 5;
      }
      // counts
      map<long, short int> counts;
      map<long, vector<long>> tmpMap;
      tmpMap = this->data[chromosome][value];
      // init data counts
      for (long j = start; j <= end; j++) {
        counts[j] = 0;
      }
      // get pair of iterators on range of intervals including the start point
      auto lb = tmpMap.lower_bound(start);
      auto ub = tmpMap.upper_bound(end);
      // cerr << "\t\t\tGet range of starts from : \n";
      map<long, vector<long>>::iterator myIterC;
      for (myIterC = lb; myIterC != ub; myIterC++) {
        vector<long>::iterator myIterD;
        for (myIterD = (myIterC->second).begin();
              myIterD != (myIterC->second).end(); myIterD++) {
          // cerr <<  "\t\t\t\t" << myIterC->first << " to " << *myIterD <<
          // endl;
          long i = myIterC->first;
          while (i <= *myIterD) {
            if (i > end) {
              break;
            }
            counts[i]++;
            i++;
          }
        }
      }
      vector<double> m;
      for (auto &s : counts) {
        m.push_back((double)(s.second));
      }
      double mean = moyenne_calculator(m);

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
            while (getline(issInfo, infomot , ';')) {
              if (infomot.find("SVTYPE=") == 0) {
                svtype = parseOnSep(infomot, "=")[1];
                continue;
              }
              if (infomot.find("CIEND=") == 0) {
                ciend = parseOnSep(infomot, "=")[1];
                continue;
              }
              if (infomot.find("VALUE=") == 0) {
                value = parseOnSep(infomot, "=")[1];
                continue;
              }
              outStream << infomot << ";";
            }
            outStream << "CIEND=" << ciend << ";VALUE=" << value << ";SVTYPE="; 
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
}

short int cnvCompare::getNbFile() { return this->nbFile; }

string cnvCompare::getInputFile() {
  return this->inputFile; 
}

string cnvCompare::getControlFile() {
  return this->controlFile;
}


int cnvCompare::fillMap(string incFile, string status) {
  ifstream allStream(incFile.c_str());
  string ligne;
  if (allStream) {
    while (getline(allStream, ligne)) {
      ifstream cnvStream(ligne.c_str());
      cerr << "\tadding file " << ligne << " as " << status << endl;
      if (!IsFileReadable(ligne)) {
        cerr << "File provided as " << status << " : " << ligne << " is not accessible : passsing file" << endl;
        continue;
      }

      this->nbFile++;
      this->fileMap[ligne] = status;
    }
  }
  return this->nbFile;
}


int cnvCompare::getFilterSize() {
  return this->filterSize;
}

//

void cnvCompare::setFormat(bool incVCFChoice, bool incBEDChoice) {
  this->useVCFFormat = incVCFChoice;
  this->useBEDFormat = incBEDChoice;
}

string cnvCompare::getFormat() {
  if (this->useVCFFormat) {
    return "VCF"; 
  } else {
    return "BED";
  }


}