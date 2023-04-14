all:
	g++ -c cnvCompare.cpp -O3 -std=c++2a
	g++ -c utils.cpp -O3
	g++ -c main.cpp -O3
	g++ main.o utils.o cnvCompare.o -o cnvCompare -lboost_program_options -O3

ccub:
	g++ -c cnvCompare.cpp -O3 -std=c++2a -I /work/gad/shared/bin/gcc/gcc-12.2.0/libstdc++-v3/include
	g++ -c utils.cpp -O3
	g++ -c main.cpp -O3
	g++ main.o utils.o cnvCompare.o -o cnvCompare -lboost_program_options -L /work/gad/shared/bin/lib/lib64/lib/ -O3

clean:
	rm *.o
	rm cnvCompare
