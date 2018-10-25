
all: main gen print eql


main: main.cpp 
	g++ -o main main.cpp /usr/lib/x86_64-linux-gnu/libpapi.so.5

gen: gen.cpp 
	g++ -o gen gen.cpp

print: print.cpp
	g++ -o print print.cpp

eql: eql.cpp
	g++ -o eql eql.cpp



test: 
	./main testA.dat testB.dat testC.dat 0 0
	./eql testC.dat testREZ.dat
	./main testA.dat testB.dat testC.dat 0 1
	./eql testC.dat testREZ.dat
	./main testA.dat testB.dat testC.dat 1 0
	./eql testC.dat testREZ.dat
	./main testA.dat testB.dat testC.dat 1 1
	./eql testC.dat testREZ.dat
	
proverka:
	./gen f 2000 2000 A.dat
	./gen f 2000 2000 B.dat
	./main A.dat B.dat C.dat 0 0
	./main A.dat B.dat C.dat 0 1
	./main A.dat B.dat C.dat 1 0
	./main A.dat B.dat C.dat 1 1

	


clean:
	rm -f *.o main gen print eql \
		A.dat B.dat C.dat testC.dat
	rm -f *.o commands
	rm -f *.txt
	rm -f *.png

com:
	./main A.dat B.dat C.dat 0 0
	./main A.dat B.dat C.dat 1 0
	./main A.dat B.dat C.dat 1 1

report:
	./gen f 1000 1000 A.dat
	./gen f 1000 1000 B.dat
	make com
	./gen f 2000 2000 A.dat
	./gen f 2000 2000 B.dat
	make com
	./gen f 3000 3000 A.dat
	./gen f 3000 3000 B.dat
	make com
	./gen f 4000 4000 A.dat
	./gen f 4000 4000 B.dat
	make com
	./gen f 5000 5000 A.dat
	./gen f 5000 5000 B.dat
	make com



		echo "set terminal png size 1024, 768 \n set output \"time(size).png\" \n set xlabel \"Size of matrix\" \n set ylabel \"Time (sec)\" \n set grid \n plot \"DATA32ijk.txt\" using 1:4 with lines lt rgb \"red\" lw 2 title \"Size of block 32 ijk\" , \"DATA32ikj.txt\" using 1:4 with lines lt rgb \"blue\" lw 2 title \"Size of block 32 ikj\" ,  \"DATA50ikj.txt\" using 1:4 with lines lt rgb \"green\" lw 2 title \"Size of block 50 ikj\" \nexit" > commands; gnuplot commands

		echo "set terminal png size 1024, 768 \n set output \"L1misses(size).png\" \n set xlabel \"Size of matrix\" \n set ylabel \"L1misses\" \n set grid \n plot \"DATA32ijk.txt\" using 1:5 with lines lt rgb \"red\" lw 2 title \"Size of block 32 ijk\" , \"DATA32ikj.txt\" using 1:5 with lines lt rgb \"blue\" lw 2 title \"Size of block 32 ikj\" ,  \"DATA50ikj.txt\" using 1:5 with lines lt rgb \"green\" lw 2 title \"Size of block 50 ikj\" \nexit" > commands; gnuplot commands

		echo "set terminal png size 1024, 768 \n set output \"L2misses(size).png\" \n set xlabel \"Size of matrix\" \n set ylabel \"L2misses\" \n set grid \n plot \"DATA32ijk.txt\" using 1:6 with lines lt rgb \"red\" lw 2 title \"Size of block 32 ijk\" , \"DATA32ikj.txt\" using 1:6 with lines lt rgb \"blue\" lw 2 title \"Size of block 32 ikj\" ,  \"DATA50ikj.txt\" using 1:6 with lines lt rgb \"green\" lw 2 title \"Size of block 50 ikj\" \nexit" > commands; gnuplot commands

		echo "set terminal png size 1024, 768 \n set output \"CPUcycles(size).png\" \n set xlabel \"Size of matrix\" \n set ylabel \"CPUcycles\" \n set grid \n plot \"DATA32ijk.txt\" using 1:7 with lines lt rgb \"red\" lw 2 title \"Size of block 32 ijk\" , \"DATA32ikj.txt\" using 1:7 with lines lt rgb \"blue\" lw 2 title \"Size of block 32 ikj\" ,  \"DATA50ikj.txt\" using 1:7 with lines lt rgb \"green\" lw 2 title \"Size of block 50 ikj\" \nexit" > commands; gnuplot commands

	rm -f *.o commands
	rm -f *.txt

