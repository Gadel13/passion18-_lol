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

	(echo "H 2" 1>input);
	(mpirun --oversubscribe -n 4 ./Task4 5 0 test_in.dat REZ_H_2.dat < input) 1>outfiles;
	(./eql true_H_2.dat REZ_H_2.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -e "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ* outfiles

testnH:
	make clean;
	make;
	cp ./test/*.dat ./;

	(echo "nH" 1>input);
	(mpirun --oversubscribe -n 4 ./Task4 5 0 test_in.dat REZ_nH.dat < input) 1>outfiles;
	(./eql true_nH.dat REZ_nH.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -e "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ* outfiles

testNOT:
	make clean;
	make;
	cp ./test/*.dat ./

	(echo "NOT 2" 1>input);
	(mpirun --oversubscribe  -n 4 ./Task4 5 0 test_in.dat REZ_NOT_2.dat < input) 1>outfiles;
	(./eql true_NOT_2.dat REZ_NOT_2.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -e "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ* outfiles

testCNOT:
	make clean;
	make;
	cp ./test/*.dat ./;

	(echo "CNOT 1 4" 1>input);
	(mpirun --oversubscribe  -n 4 ./Task4 5 0 test_in.dat REZ_CNOT_1_4.dat < input) 1>outfiles;
	(./eql true_CNOT_1_4.dat REZ_CNOT_1_4.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -e "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ* outfiles

testRw:
	make clean;
	make;
	cp ./test/*.dat ./;

	(echo "Rw 3 0.8" 1>input);
	(mpirun --oversubscribe -n 4 ./Task4 5 0 test_in.dat REZ_Rw_3_08.dat < input) 1>outfiles;
	(./eql true_Rw_3_08.dat REZ_Rw_3_08.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -e "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ* outfiles

testCRw:
	make clean;
	make;
	cp ./test/*.dat ./;

	(echo "CRw 2 3 0.5" 1>input);
	(mpirun --oversubscribe  -n 4 ./Task4 5 0 test_in.dat REZ_CRw_2_3_05.dat < input) 1>outfiles;
	(./eql true_CRw_2_3_05.dat REZ_CRw_2_3_05.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -e "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ* outfiles

cpplint:
	cpplint --linelength=150 --filter=-legal/copyright,-build/include_subdir,-build/namespaces,-readability/casting,-runtime,-build/include_order Task4.cpp && echo -e "\033[1;36m Task4.cpp cpplint done\033[0m" || echo -e "\033[37;1;41m Task4.cpp cpplint fail \033[0m";
	cpplint --linelength=150 --filter=-legal/copyright,-build/include_subdir,-build/namespaces,-readability/casting,-runtime,-build/header_guard,-build/include_order gates.h && echo -e "\033[1;36m gates.h cpplint done\033[0m" || echo -e "\033[37;1;41m gates.h cpplint fail \033[0m";
	cpplint --linelength=150 --filter=-legal/copyright,-build/include_subdir,-build/namespaces,-readability/casting,-runtime,-build/include_order gates.cpp && echo -e "\033[1;36m gates.cpp cpplint done\033[0m" || echo -e "\033[37;1;41m gates.cpp cpplint fail \033[0m";
	cpplint --linelength=150 --filter=-legal/copyright,-build/include_subdir,-build/namespaces,-readability/casting,-runtime,-build/include_order eql.cpp && echo -e "\033[1;36m eql.cpp cpplint done\033[0m" || echo -e "\033[37;1;41m eql.cpp cpplint fail \033[0m";





# without --oversubscribe
testm:
	make testHm
	make testnHm
	make testNOTm
	make testCNOTm
	make testRwm
	make testCRwm

testHm:
	make clean;
	make;
	cp ./test/*.dat ./;

	(echo "H 2" 1>input);
	(mpirun -n 4 ./Task4 5 0 test_in.dat REZ_H_2.dat < input) 1>outfiles;
	(./eql true_H_2.dat REZ_H_2.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -e "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ* outfiles

testnHm:
	make clean;
	make;
	cp ./test/*.dat ./;

	(echo "nH" 1>input);
	(mpirun -n 4 ./Task4 5 0 test_in.dat REZ_nH.dat < input) 1>outfiles;
	(./eql true_nH.dat REZ_nH.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -e "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ* outfiles

testNOTm:
	make clean;
	make;
	cp ./test/*.dat ./

	(echo "NOT 2" 1>input);
	(mpirun -n 4 ./Task4 5 0 test_in.dat REZ_NOT_2.dat < input) 1>outfiles;
	(./eql true_NOT_2.dat REZ_NOT_2.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -e "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ* outfiles

testCNOTm:
	make clean;
	make;
	cp ./test/*.dat ./;

	(echo "CNOT 1 4" 1>input);
	(mpirun -n 4 ./Task4 5 0 test_in.dat REZ_CNOT_1_4.dat < input) 1>outfiles;
	(./eql true_CNOT_1_4.dat REZ_CNOT_1_4.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -e "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ* outfiles

testRwm:
	make clean;
	make;
	cp ./test/*.dat ./;

	(echo "Rw 3 0.8" 1>input);
	(mpirun -n 4 ./Task4 5 0 test_in.dat REZ_Rw_3_08.dat < input) 1>outfiles;
	(./eql true_Rw_3_08.dat REZ_Rw_3_08.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -e "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ* outfiles

testCRwm:
	make clean;
	make;
	cp ./test/*.dat ./;

	(echo "CRw 2 3 0.5" 1>input);
	(mpirun -n 4 ./Task4 5 0 test_in.dat REZ_CRw_2_3_05.dat < input) 1>outfiles;
	(./eql true_CRw_2_3_05.dat REZ_CRw_2_3_05.dat && echo -e "\033[1;36mEQ files!\033[0m" || echo -e "\033[37;1;41m NOT_EQ_files! \033[0m");
	rm input;
	rm true* REZ* outfiles

