
all: main 


main: main.cpp 
	mpic++ -std=c++11 -O3 -o main ./main.cpp


clean:
	rm -f *.o main
	rm -f commands
	rm -f *.txt
	rm -f *.png


report:
	mpirun -n 2 ./main 2 100000000 rez.txt
	mpirun -n 4 ./main 2 100000000 rez.txt
	mpirun -n 6 ./main 2 100000000 rez.txt
	mpirun -n 8 ./main 2 100000000 rez.txt
	mpirun -n 12 ./main 2 100000000 rez.txt
	mpirun -n 16 ./main 2 100000000 rez.txt
	mpirun -n 24 ./main 2 100000000 rez.txt
	mpirun -n 32 ./main 2 100000000 rez.txt
	mpirun -n 64 ./main 2 100000000 rez.txt






		echo "set terminal png size 1024, 768 \n set output \"totaltime(num proc).png\" \n set xlabel \"Process num\" \n set ylabel \"Total Time (sec)\" \n set grid \n plot \"TT.txt\" using 1:2 with lines lt rgb \"red\" lw 2 title \"Total time(cpu num) 2-100000000\"  \nexit" > commands; gnuplot commands

		echo "set terminal png size 1024, 768 \n set output \"maxtime(num proc).png\" \n set xlabel \"Process num\" \n set ylabel \"Max Time (sec)\" \n set grid \n plot \"MT.txt\" using 1:2 with lines lt rgb \"red\" lw 2 title \"Max time(cpu num) 2-100000000\"  \nexit" > commands; gnuplot commands



	rm -f *.o commands
	rm -f *.txt

