
all: main gen print speedup


main: main.cpp 
	mpic++ -std=c++11 -O3 -o main ./main.cpp

gen: gen.cpp 
	g++ -std=c++11 -O3 -o gen gen.cpp

print: print.cpp
	g++ -std=c++11 -O3 -o print print.cpp

speedup: speedup.cpp
	g++ -std=c++11 -O3 -o speedup speedup.cpp


clean:
	rm -f *.o main gen print
	rm -f *.dat
	rm -f *.o commands
	rm -f *.o *.png *.txt


test0:
	make clean
	make
	./gen d 500 300 A.dat
	./gen d 300 1 B.dat
	mpirun -n 8 ./main A.dat B.dat C.dat

test1:
	make clean
	make
	./gen d 2 3 A.dat
	./gen d 3 1 B.dat
	mpirun -n 8 ./main A.dat B.dat C.dat

report:
	make clean
	make
	./gen d 2048 2048 A.dat
	./gen d 2048 1 B.dat
	mpirun -n 1 ./main A.dat B.dat C.dat
	./gen d 2048 2048 A.dat
	./gen d 2048 1 B.dat
	mpirun -n 2 ./main A.dat B.dat C.dat
	./gen d 2048 2048 A.dat
	./gen d 2048 1 B.dat
	mpirun -n 4 ./main A.dat B.dat C.dat
	./gen d 2048 2048 A.dat
	./gen d 2048 1 B.dat
	mpirun -n 8 ./main A.dat B.dat C.dat
	./gen d 2048 2048 A.dat
	./gen d 2048 1 B.dat
	mpirun -n 16 ./main A.dat B.dat C.dat
	./gen d 2048 2048 A.dat
	./gen d 2048 1 B.dat
	mpirun -n 32 ./main A.dat B.dat C.dat

	./speedup MT.txt

	echo "set terminal png size 1024, 768 \n set output \"totaltime(num proc).png\" \n set xlabel \"Process num\" \n set ylabel \"Total Time (sec)\" \n set grid \n plot \"TT.txt\" using 1:2 with lines lt rgb \"red\" lw 2 title \"Total time(cpu num)\"  \nexit" > commands; gnuplot commands

	echo "set terminal png size 1024, 768 \n set output \"maxtime(num proc).png\" \n set xlabel \"Process num\" \n set ylabel \"Max Time (sec)\" \n set grid \n plot \"MT.txt\" using 1:2 with lines lt rgb \"red\" lw 2 title \"Max time(cpu num)\"  \nexit" > commands; gnuplot commands



	rm -f commands
	
