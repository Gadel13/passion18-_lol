PROG = Task5
CFLAGS = -Werror -Wall -fopenmp -fstack-protector

.PHONY: test

all: gates.o Task5 eql conv_to_binary

eql: eql.cpp
	g++ -Werror -Wall -o eql eql.cpp

conv_to_binary: conv_to_binary.cpp
	g++ -Werror -Wall -o conv_to_binary conv_to_binary.cpp

gates.o: gates.cpp gates.h 
	mpicxx $(CFLAGS) -c  $<  -o  $@

$(PROG): Task5.cpp gates.o
	mpicxx $(CFLAGS) Task5.cpp gates.o -o $(PROG)


clean:
	rm -f *.o F Task5 eql conv_to_binary *.dat *.txt input test_* REZ*

test:
	make test1
	make test2
	make test4
	make test8
	make clean

test1:
	make clean;
	make;
	cp ./test/* ./;
	./conv_to_binary test_in test_in.dat
	./conv_to_binary test_true test_true.dat

	mpirun --oversubscribe -n 1 ./Task5 3 0 test_in.dat REZ.dat;
	(./eql test_true.dat REZ.dat && echo -e "\033[1;36m n1 TEST QFT EQ files!\033[0m") || (echo -e "\033[37;1;41m n1 TEST QFT NOT_EQ_files! \033[0m" && exit 1);
	rm test_* REZ*

test2:
	make clean;
	make;
	cp ./test/* ./;
	./conv_to_binary test_in test_in.dat
	./conv_to_binary test_true test_true.dat

	mpirun --oversubscribe -n 2 ./Task5 3 0 test_in.dat REZ.dat;
	(./eql test_true.dat REZ.dat && echo -e "\033[1;36m n2 TEST QFT EQ files!\033[0m") || (echo -e "\033[37;1;41m n2 TEST QFT NOT_EQ_files! \033[0m" && exit 1);
	rm test_* REZ*

test4:
	make clean;
	make;
	cp ./test/* ./;
	./conv_to_binary test_in test_in.dat
	./conv_to_binary test_true test_true.dat

	mpirun --oversubscribe -n 4 ./Task5 3 0 test_in.dat REZ.dat;
	(./eql test_true.dat REZ.dat && echo -e "\033[1;36m n4 TEST QFT EQ files!\033[0m") || (echo -e "\033[37;1;41m n4 TEST QFT NOT_EQ_files! \033[0m" && exit 1);
	rm test_* REZ*

test8:
	make clean;
	make;
	cp ./test/* ./;
	./conv_to_binary test_in test_in.dat
	./conv_to_binary test_true test_true.dat

	mpirun --oversubscribe -n 8 ./Task5 3 0 test_in.dat REZ.dat;
	(./eql test_true.dat REZ.dat && echo -e "\033[1;36m n8 TEST QFT EQ files!\033[0m") || (echo -e "\033[37;1;41m n8 TEST QFT NOT_EQ_files! \033[0m" && exit 1);
	rm test_* REZ*

cpplint:
	(cpplint --linelength=150 --filter=-legal/copyright,-build/include_subdir,-build/namespaces,-readability/casting,-runtime,-build/include_order,-build/include Task5.cpp && echo -e "\033[1;36m Task5.cpp cpplint done\033[0m") || (echo -e "\033[37;1;41m Task5.cpp cpplint fail \033[0m" && exit 1);
	(cpplint --linelength=150 --filter=-legal/copyright,-build/include_subdir,-build/namespaces,-readability/casting,-runtime,-build/header_guard,-build/include_order,-build/include gates.h && echo -e "\033[1;36m gates.h cpplint done\033[0m") || (echo -e "\033[37;1;41m gates.h cpplint fail \033[0m" && exit 1);
	(cpplint --linelength=150 --filter=-legal/copyright,-build/include_subdir,-build/namespaces,-readability/casting,-runtime,-build/include_order,-build/include gates.cpp && echo -e "\033[1;36m gates.cpp cpplint done\033[0m") || (echo -e "\033[37;1;41m gates.cpp cpplint fail \033[0m" && exit 1);
	(cpplint --linelength=150 --filter=-legal/copyright,-build/include_subdir,-build/namespaces,-readability/casting,-runtime,-build/include_order,-build/include eql.cpp && echo -e "\033[1;36m eql.cpp cpplint done\033[0m") || (echo -e "\033[37;1;41m eql.cpp cpplint fail \033[0m" && exit 1);
	(cpplint --linelength=150 --filter=-legal/copyright,-build/include_subdir,-build/namespaces,-readability/casting,-runtime,-build/include_order,-build/include conv_to_binary.cpp && echo -e "\033[1;36m conv_to_binary.cpp cpplint done\033[0m") || (echo -e "\033[37;1;41m conv_to_binary.cpp cpplint fail \033[0m" && exit 1);





# without --oversubscribe
testm:
	make test1m
	make test2m
	make test4m
	make test8m
	make clean

test1m:
	make clean;
	make;
	cp ./test/* ./;
	./conv_to_binary test_in test_in.dat
	./conv_to_binary test_true test_true.dat

	mpirun -n 1 ./Task5 3 0 test_in.dat REZ.dat;
	(./eql test_true.dat REZ.dat && echo -e "\033[1;36m TEST n1 QFT EQ files!\033[0m") || (echo -e "\033[37;1;41m TEST n1 QFT NOT_EQ_files! \033[0m" && exit 1);
	rm test_* REZ*

test2m:
	make clean;
	make;
	cp ./test/* ./;
	./conv_to_binary test_in test_in.dat
	./conv_to_binary test_true test_true.dat

	mpirun -n 2 ./Task5 3 0 test_in.dat REZ.dat;
	(./eql test_true.dat REZ.dat && echo -e "\033[1;36m TEST n2 QFT EQ files!\033[0m") || (echo -e "\033[37;1;41m TEST n2 QFT NOT_EQ_files! \033[0m" && exit 1);
	rm test_* REZ*

test4m:
	make clean;
	make;
	cp ./test/* ./;
	./conv_to_binary test_in test_in.dat
	./conv_to_binary test_true test_true.dat

	mpirun -n 4 ./Task5 3 0 test_in.dat REZ.dat;
	(./eql test_true.dat REZ.dat && echo -e "\033[1;36m TEST n4 QFT EQ files!\033[0m") || (echo -e "\033[37;1;41m TEST n4 QFT NOT_EQ_files! \033[0m" && exit 1);
	rm test_* REZ*

test8m:
	make clean;
	make;
	cp ./test/* ./;
	./conv_to_binary test_in test_in.dat
	./conv_to_binary test_true test_true.dat

	mpirun -n 8 ./Task5 3 0 test_in.dat REZ.dat;
	(./eql test_true.dat REZ.dat && echo -e "\033[1;36m TEST n8 QFT EQ files!\033[0m") || (echo -e "\033[37;1;41m TEST n8 QFT NOT_EQ_files! \033[0m" && exit 1);
	rm test_* REZ*

