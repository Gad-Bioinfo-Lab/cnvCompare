#ifndef _CNVCOMPARE_H
#define _CNVCOMPARE_H

#include <vector>
#include <map>

class cnvCompare
{
	public:
		cnvCompare();
		cnvCompare( std::string , int );
		cnvCompare(std::string , std::string, int);
		void mainLoop();
		void altLoop();
		void getData();
		void cleanData();
		void getDataWhole(std::string);
		void computeCounts();
		void computeChrCounts(std::string);
		short int getNbFile();

	private:
		std::string inputFile;
		std::string controlFile;
		int nbThread;
		std::map< std::string , std::map< unsigned int , std::map< long , std::vector< long > > > > data;
		std::map<unsigned int , std::map<long, short> > dataByChr;
		std::vector<std::string> chromosomeMap ;
		short int nbFile = 0;
		bool useControls = false; 

};
#endif
