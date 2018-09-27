
all: main gen print eql


main: main.cpp 
	g++ -o main main.cpp

gen: gen.cpp 
	g++ -o gen gen.cpp

print: print.cpp
	g++ -o print print.cpp

eql: eql.cpp
	g++ -o eql eql.cpp

test: 
	./main testA.dat testB.dat testC0.dat 0
	./eql testC0.dat testREZ.dat
	./main testA.dat testB.dat testC1.dat 1
	./eql testC1.dat testREZ.dat
	./main testA.dat testB.dat testC2.dat 2
	./eql testC2.dat testREZ.dat
	./main testA.dat testB.dat testC3.dat 3
	./eql testC3.dat testREZ.dat
	./main testA.dat testB.dat testC4.dat 4
	./eql testC4.dat testREZ.dat
	./main testA.dat testB.dat testC5.dat 5
	./eql testC5.dat testREZ.dat
	


clean:
	rm -f *.o main
	rm -f *.o gen
	rm -f *.o print
	rm -f *.o eql
	rm -f *.o A.dat
	rm -f *.o B.dat
	rm -f *.o C.dat
	rm -f *.o commands
	rm -f *.o DATA.txt
	rm -f *.o testC0.dat
	rm -f *.o testC1.dat
	rm -f *.o testC2.dat
	rm -f *.o testC3.dat
	rm -f *.o testC4.dat
	rm -f *.o testC5.dat
	rm -f *.o graph50x50.png
	rm -f *.o graph200x200.png
	rm -f *.o graph300x300.png
	rm -f *.o graph500x500.png
	rm -f *.o graph1000x1000.png


report:
	rm -f *.o DATA.txt
	./gen d 50 50 A.dat
	./gen d 50 50 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	 echo "set terminal png size 1024, 768 \nset output 'graph50x50.png' \nset xtics (\"ijk\" 0, \"ikj\" 1, \"kij\" 2, \"jik\" 3, \"jki\" 4, \"kji\" 5) \nplot 'DATA.txt' u 1:2 w linesp \nexit" > commands; gnuplot commands

	rm -f *.o DATA.txt
	./gen d 200 200 A.dat
	./gen d 200 200 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	echo "set terminal png size 1024, 768 \nset output 'graph200x200.png' \nset xtics (\"ijk\" 0, \"ikj\" 1, \"kij\" 2, \"jik\" 3, \"jki\" 4, \"kji\" 5) \nplot 'DATA.txt' u 1:2 w linesp \nexit" > commands; gnuplot commands
	
	rm -f *.o DATA.txt
	./gen d 300 300 A.dat
	./gen d 300 300 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	echo "set terminal png size 1024, 768 \nset output 'graph300x300.png' \nset xtics (\"ijk\" 0, \"ikj\" 1, \"kij\" 2, \"jik\" 3, \"jki\" 4, \"kji\" 5) \nplot 'DATA.txt' u 1:2 w linesp \nexit" > commands; gnuplot commands

	rm -f *.o DATA.txt
	./gen d 500 500 A.dat
	./gen d 500 500 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	echo "set terminal png size 1024, 768 \nset output 'graph500x500.png' \nset xtics (\"ijk\" 0, \"ikj\" 1, \"kij\" 2, \"jik\" 3, \"jki\" 4, \"kji\" 5) \nplot 'DATA.txt' u 1:2 w linesp \nexit" > commands; gnuplot commands

	rm -f *.o DATA.txt
	./gen d 1000 1000 A.dat
	./gen d 1000 1000 B.dat
	./main A.dat B.dat C.dat 0
	./main A.dat B.dat C.dat 1
	./main A.dat B.dat C.dat 2
	./main A.dat B.dat C.dat 3
	./main A.dat B.dat C.dat 4
	./main A.dat B.dat C.dat 5
	echo "set terminal png size 1024, 768 \nset output 'graph1000x1000.png' \nset xtics (\"ijk\" 0, \"ikj\" 1, \"kij\" 2, \"jik\" 3, \"jki\" 4, \"kji\" 5) \nplot 'DATA.txt' u 1:2 w linesp \nexit" > commands; gnuplot commands
	
