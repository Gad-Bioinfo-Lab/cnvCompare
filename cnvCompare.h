#ifndef _CNVCOMPARE_H
#define _CNVCOMPARE_H

#include <vector>
#include <map>

class cnvCompare
{
	public:
		cnvCompare();
		cnvCompare( std::string , int , int);
		cnvCompare(std::string , std::string, int, int);
		void mainLoop();
		void altLoop();
		void getData();
		void cleanData();
		void getDataWhole(std::string);
		void computeCounts();
		void computeChrCounts(std::string);
		short int getNbFile();
		int fillMap(std::string, std::string); 
		std::string getControlFile();
		std::string getInputFile();
		std::string getFormat(); 
		int getFilterSize();
		void setFormat(bool, bool);
		std::vector<std::string> parseBEDLine(std::string);
		std::vector<std::string> parseVCFLine(std::string);

	private:
		std::string inputFile;
		std::string controlFile;
		int nbThread;
		std::map< std::string , std::map< unsigned int , std::map< long , std::vector< long > > > > data;
		std::map<unsigned int , std::map<long, short> > dataByChr;
		std::vector<std::string> chromosomeMap ;
		short int nbFile = 0;
		bool useControls = false; 
		std::map<std::string, std::string> fileMap; 
		int filterSize;
		bool useVCFFormat; 
		bool useBEDFormat;

};
#endif
