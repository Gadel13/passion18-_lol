
all: main gen print speedup random


main: main.cpp 
	mpicxx  -O3 -o main ./main.cpp

gen: gen.cpp 
	g++ -O3 -o gen gen.cpp

print: print.cpp
	g++  -O3 -o print print.cpp

speedup: speedup.cpp
	g++ -O3 -o speedup speedup.cpp

random: random.cpp
	g++ -o random random.cpp



clean:
	rm -f *.o main gen print speedup random eql
	rm -f A.dat B.dat C.dat testC.dat
	rm -f *.o commands
	rm -f *.o *.png *.txt


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
	
