all:
	g++ -c cnvCompare.cpp -O3
	g++ -c utils.cpp -O3
	g++ -c main.cpp -O3
	g++ main.o utils.o cnvCompare.o -o cnvCompare -lboost_program_options

ccub:
	g++ -c cnvCompare.cpp -O3 -std=gnu++0x
	g++ -c utils.cpp -O3
	g++ -c main.cpp -O3
	g++ main.o utils.o cnvCompare.o -o cnvCompare -lhts -lboost_program_options 

clean:
	rm *.o
	rm cnvCompare
