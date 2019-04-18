PROG = Task4
CFLAGS = -Werror -Wall -fopenmp

.PHONY: test

all: gates.o Task4 

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
	make clean
	make
	cp ./test/*.dat ./

	echo "H\n5" 1>input
	mpirun --oversubscribe -n 4 ./Task4 10 0 test_in.dat REZ_H_5.dat < input
	cmp true_H_5.dat REZ_H_5.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -en "\033[37;1;41m NOT_EQ_files! \033[0m"
	rm input
	rm true* REZ*

testnH:
	make clean
	make
	cp ./test/*.dat ./

	echo "nH" 1>input
	mpirun --oversubscribe  -n 4 ./Task4 10 0 test_in.dat REZ_nH.dat < input
	cmp true_nH.dat REZ_nH.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -en "\033[37;1;41m NOT_EQ_files! \033[0m"
	rm input
	rm true* REZ*

testNOT:
	make clean
	make
	cp ./test/*.dat ./

	echo "NOT\n5" 1>input
	mpirun --oversubscribe -n 4 ./Task4 10 0 test_in.dat REZ_NOT_5.dat < input
	cmp true_NOT_5.dat REZ_NOT_5.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -en "\033[37;1;41m NOT_EQ_files! \033[0m"
	rm input
	rm true* REZ*

testCNOT:
	make clean
	make
	cp ./test/*.dat ./

	echo "CNOT\n2\n5" 1>input
	mpirun  --oversubscribe -n 4 ./Task4 10 0 test_in.dat REZ_CNOT_2_5.dat < input
	cmp true_CNOT_2_5.dat REZ_CNOT_2_5.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -en "\033[37;1;41m NOT_EQ_files! \033[0m"
	rm input
	rm true* REZ*

testRw:
	make clean
	make
	cp ./test/*.dat ./

	echo "Rw\n5\n0.5" 1>input
	mpirun --oversubscribe -n 4 ./Task4 10 0 test_in.dat REZ_Rw_5_05.dat < input
	cmp true_Rw_5_05.dat REZ_Rw_5_05.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -en "\033[37;1;41m NOT_EQ_files! \033[0m"
	rm input
	rm true* REZ*

testCRw:
	make clean
	make
	cp ./test/*.dat ./

	echo "CRw\n2\n5\n0.5" 1>input
	mpirun --oversubscribe -n 4 ./Task4 10 0 test_in.dat REZ_CRw_2_5_05.dat < input
	cmp true_CRw_2_5_05.dat REZ_CRw_2_5_05.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -en "\033[37;1;41m NOT_EQ_files! \033[0m"
	rm input

	rm true* REZ*



cpplint:
	cpplint --linelength=125 --filter=-legal/copyright,-readability/braces,-build/include_subdir,-build/namespaces,-readability/casting,-whitespace/braces,-runtime Task4.cpp
