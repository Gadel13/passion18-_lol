
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
	./main testA.txt testB.txt testC0.txt 0
	./eql testC0.txt testREZ.txt
	./main testA.txt testB.txt testC1.txt 1
	./eql testC1.txt testREZ.txt
	./main testA.txt testB.txt testC2.txt 2
	./eql testC2.txt testREZ.txt
	./main testA.txt testB.txt testC3.txt 3
	./eql testC3.txt testREZ.txt
	./main testA.txt testB.txt testC4.txt 4
	./eql testC4.txt testREZ.txt
	./main testA.txt testB.txt testC5.txt 5
	./eql testC5.txt testREZ.txt
	


clean:
	rm -f *.o main
	rm -f *.o gen
	rm -f *.o print

report:
	gnuplot commands
	
