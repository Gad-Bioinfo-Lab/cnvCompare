all:
	g++ -c cnvCompare.cpp -O3
	g++ -c utils.cpp -O3
	g++ -c main.cpp -O3
	g++ main.o utils.o cnvCompare.o -o cnvCompare -lboost_program_options

ccub:
	g++ -c cnvCompare.cpp -O3 -std=gnu++0x -I /usr/local/src/htslib-develop -I /work/gad/shared/bin/htslib-1.9
	g++ -c utils.cpp -O3
	g++ -c main.cpp -O3 -I /usr/local/src/htslib-develop -I /work/gad/shared/bin/htslib-1.9
	g++ main.o utils.o cnvCompare.o -o cnvCompare -lhts -lboost_program_options -L /work/gad/shared/bin/htslib-1.9

clean:
	rm *.o
	rm cnvCompare
