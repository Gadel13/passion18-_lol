
all: Task3 F avr


Task3: Task3.cpp 
	mpicxx   -fopenmp -o Task3 ./Task3.cpp

F: F.cpp
	g++ -o F F.cpp

avr: avr.cpp
	g++ -o avr avr.cpp

clean:
	rm -f *.o F Task3 avr data.txt *.dat *.txt

test23:
	mpirun -n 8 ./Task3 23 1 1 n23.dat
	mpirun -n 8 ./Task3 23 0 0 n23.dat
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt

	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt

	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt

	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt

	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt

	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt

	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt

	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt

	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt
	mpirun -n 8 ./Task3 23 0 1 n23.dat
	./F 23 REZ.dat REZideal.dat F23.txt

test24:
	mpirun -n 8 ./Task3 24 1 1 n24.dat
	mpirun -n 8 ./Task3 24 0 0 n24.dat
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt

	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt

	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt

	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt

	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt

	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt

	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt

	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt

	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt
	mpirun -n 8 ./Task3 24 0 1 n24.dat
	./F 24 REZ.dat REZideal.dat F24.txt

test25:
	mpirun -n 8 ./Task3 25 1 1 n25.dat
	mpirun -n 8 ./Task3 25 0 0 n25.dat
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt

	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt

	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt

	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt

	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt

	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt

	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt

	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt

	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt
	mpirun -n 8 ./Task3 25 0 1 n25.dat
	./F 25 REZ.dat REZideal.dat F25.txt

test26:
	mpirun -n 8 ./Task3 26 1 1 n26.dat
	mpirun -n 8 ./Task3 26 0 0 n26.dat
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt

	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt

	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt

	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt

	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt

	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt

	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt

	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt

	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt
	mpirun -n 8 ./Task3 26 0 1 n26.dat
	./F 26 REZ.dat REZideal.dat F26.txt


test27:
	mpirun -n 8 ./Task3 27 1 1 n27.dat
	mpirun -n 8 ./Task3 27 0 0 n27.dat
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt

	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt

	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt

	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt

	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt

	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt

	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt

	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt

	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
	mpirun -n 8 ./Task3 27 0 1 n27.dat
	./F 27 REZ.dat REZideal.dat F27.txt
