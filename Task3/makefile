
all: Task3 avr


Task3: Task3.cpp 
	mpicxx   -fopenmp -o Task3 ./Task3.cpp


avr: avr.cpp
	g++ -o avr avr.cpp

clean:
	rm -f *.o F Task3 avr data.txt *.dat *.txt


