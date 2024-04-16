EXE = cnvCompare
OBJ = cnvCompare.o utils.o
OBJM = main.o

# LIBS = -lboost_regex -lboost_chrono -lboost_program_options -lboost_thread -lboost_system -lpthread -lz -lhts  -lboost_log -lboost_log_setup -L /work/gad/shared/bin/htslib-1.9 -L /work/gad/shared/bin/lib -L /soft/c7/boost/1.67.0/intel/2018/lib/ -lm -limf -lirc -L /soft/c7/intel/2019/lib/intel64 
LIBS = -lboost_regex -lboost_chrono -lboost_program_options -lboost_filesystem
OPT = -O3 -Wunused-parameter -Wextra -Wunused-variable
OPTP = -Wunused-parameter -Wextra -Wunused-variable -g
INC = -std=c++20 -I src/

$(EXE):	$(OBJ) $(OBJM)
	g++ -o $(EXE) $(OBJ) $(OBJM) $(OPT) $(LIBS)

$(OBJ):%.o: %.cpp %.h
	g++ -c $< $(OPT) $(INC)

$(OBJM):%.o: %.cpp
		g++ -c $< $(OPT) $(INC)

clean:
	rm $(OBJ) $(OBJM) $(EXE)

profile:
	g++ -c src/cnvCompare.cpp $(OPTP) $(INC)
	g++ -c src/utils.cpp $(OPTP) $(INC)
	g++ -c src/main.cpp $(OPTP) $(INC)
	g++ -o cnvComparep cnvCompare.o utils.o main.o $(OPTP) $(LIBS)

local:
	g++ -c src/cnvCompare.cpp $(OPT) $(INC)
	g++ -c src/utils.cpp $(OPT) $(INC)
	g++ -c src/main.cpp $(OPT) $(INC)
	g++ -o cnvCompare cnvCompare.o utils.o main.o $(OPT) $(LIBS)
