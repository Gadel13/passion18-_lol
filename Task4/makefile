PROG = Task4
CFLAGS = -Werror -Wall -fopenmp -fstack-protector

.PHONY: test

all: gates.o Task4 eql

eql: eql.cpp
	g++ -Werror -Wall -o eql eql.cpp

gates.o: gates.cpp gates.h 
	mpicxx $(CFLAGS) -c  $<  -o  $@

$(PROG): Task4.cpp gates.o
	mpicxx $(CFLAGS) Task4.cpp gates.o -o $(PROG)


clean:
	rm -f *.o F Task4 *.dat *.txt

test:
	make testH
	make testnH
	make testNOT
	make testCNOT
	make testRw
	make testCRw

testH:
	make clean;
	make;
	cp ./test/*.dat ./;

	(echo "H\n2" 1>input);
	(mpirun  -n 4 ./Task4 5 0 test_in.dat REZ_H_2.dat < input);
	(./eql true_H_2.dat REZ_H_2.dat && echo "\033[1;36mEQ files!\033[0m" || echo  "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ*

testnH:
	make clean;
	make;
	cp ./test/*.dat ./;

	(echo "nH" 1>input);
	(mpirun  -n 4 ./Task4 5 0 test_in.dat REZ_nH.dat < input);
	(./eql true_nH.dat REZ_nH.dat && echo "\033[1;36mEQ files!\033[0m" || echo  "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ*

testNOT:
	make clean;
	make;
	cp ./test/*.dat ./

	(echo "NOT\n2" 1>input);
	(mpirun   -n 4 ./Task4 5 0 test_in.dat REZ_NOT_2.dat < input);
	(./eql true_NOT_2.dat REZ_NOT_2.dat && echo "\033[1;36mEQ files!\033[0m" || echo  "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ*

testCNOT:
	make clean;
	make;
	cp ./test/*.dat ./;

	(echo "CNOT\n1\n4" 1>input);
	(mpirun   -n 4 ./Task4 5 0 test_in.dat REZ_CNOT_1_4.dat < input);
	(./eql true_CNOT_1_4.dat REZ_CNOT_1_4.dat && echo "\033[1;36mEQ files!\033[0m" || echo  "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ*

testRw:
	make clean;
	make;
	cp ./test/*.dat ./;

	(echo "Rw\n3\n0.8" 1>input);
	(mpirun  -n 4 ./Task4 5 0 test_in.dat REZ_Rw_3_08.dat < input);
	(./eql true_Rw_3_08.dat REZ_Rw_3_08.dat && echo "\033[1;36mEQ files!\033[0m" || echo  "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ*

testCRw:
	make clean;
	make;
	cp ./test/*.dat ./;

	(echo "CRw\n2\n3\n0.5" 1>input);
	(mpirun   -n 4 ./Task4 5 0 test_in.dat REZ_CRw_2_3_05.dat < input);
	(./eql true_CRw_2_3_05.dat REZ_CRw_2_3_05.dat && echo "\033[1;36mEQ files!\033[0m" || echo  "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ*

cpplint:
	cpplint --linelength=125 --filter=-legal/copyright,-readability/braces,-build/include_subdir,-build/namespaces,-readability/casting,-whitespace/braces,-runtime Task4.cpp
