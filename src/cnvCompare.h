#ifndef _CNVCOMPARE_H
#define _CNVCOMPARE_H

#include <vector>
#include <map>
#include <unordered_map>
#include <list>

class cnvCompare
{
	public:
		cnvCompare();
		cnvCompare(std::string , int , int);
		cnvCompare(std::string , std::string, int, int);
		void mainLoop();
		void altLoop();
		void fastLoop(); 
		void getDataWhole();
		void getDatabyChr(std::string);
		void getDataFast(); 
		void computeCountsWhole();
		void computeCountsbyChr(std::string);
		void computeCountsFast();
		void cleanData();
		short int getNbFile();
		short int getNbIndividual();
		int fillMap(std::string, std::string); 
		std::string getControlFile();
		std::string getInputFile();
		std::string getFormat(); 
		int getFilterSize();
		void setFormat(bool, bool);
		std::vector<std::string> parseBEDLine(std::string);
		std::vector<std::string> parseVCFLine(std::string);
		void setSuffix(std::string);
		std::string getSuffix(); 
		std::string getDictFile(); 
		void setDictFile(std::string);
		void parseDictFile(std::string);
		void setHasDict(bool);

	private:
		std::string inputFile;
		std::string controlFile;
		int nbThread;
		std::unordered_map< std::string , std::unordered_map< unsigned int , std::unordered_map<long, short> > > data;
		std::unordered_map<unsigned int , std::unordered_map<long, short> > dataByChr;
		std::unordered_map<std::string, std::unordered_map<unsigned int, std::map<long, short> > > breakpoints; 
		std::unordered_map<std::string, std::unordered_map<unsigned int, std::map<long, short> > > invbreakpoints; 
		std::unordered_map<std::string, std::unordered_map<unsigned int, std::map<long, short> > > trnbreakpoints; 
		std::vector<std::string> chromosomeMap;
		short int nbFile = 0;
		short int nbIndividual = 0;
		bool useControls = false; 
		std::map<std::string, std::string> fileMap; 
		int filterSize;
		bool useVCFFormat; 
		bool useBEDFormat;
		int populateChr(); 
		std::string suffix;
		std::string dictFile; 
		bool hasDict; 
		void watchHeader(std::string);
};
#endif
