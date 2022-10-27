default:
	g++ -std=c++11 -fpermissive -O3 -funroll-loops -o temp main.cpp Grid.cpp VarBox.cpp Solver.cpp U.cpp V.cpp Pressure.cpp Temperature.cpp -I/home/ztdep/lib/silo-4.10.2/silo/include/ -L/home/ztdep/lib/silo-4.10.2/silo/lib64 -lsilo
