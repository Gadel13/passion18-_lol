
all: main 


main: main.cpp 
	g++ -std=c++11 -O3 -o main  main.cpp /usr/lib/x86_64-linux-gnu/libpthread.so


clean:
	rm -f *.o main
	rm -f commands
	rm -f *.txt
	rm -f *.png


report:
	make clean
	make 
	./main 2 100000000 rez.txt 2
	./main 2 100000000 rez.txt 4
	./main 2 100000000 rez.txt 6
	./main 2 100000000 rez.txt 8
	./main 2 100000000 rez.txt 12
	./main 2 100000000 rez.txt 16
	./main 2 100000000 rez.txt 24
	./main 2 100000000 rez.txt 32
	./main 2 100000000 rez.txt 64






		echo "set terminal png size 1024, 768 \n set output \"totaltime(num proc).png\" \n set xlabel \"Process num\" \n set ylabel \"Total Time (sec)\" \n set grid \n plot \"TT.txt\" using 1:2 with lines lt rgb \"red\" lw 2 title \"Total time(cpu num) 2-100000000\"  \nexit" > commands; gnuplot commands

		echo "set terminal png size 1024, 768 \n set output \"maxtime(num proc).png\" \n set xlabel \"Process num\" \n set ylabel \"Max Time (sec)\" \n set grid \n plot \"MT.txt\" using 1:2 with lines lt rgb \"red\" lw 2 title \"Max time(cpu num) 2-100000000\"  \nexit" > commands; gnuplot commands



	rm -f *.o commands
	rm -f *.txt

