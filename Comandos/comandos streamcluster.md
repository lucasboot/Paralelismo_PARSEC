#### Pthreads
g++ -pg -O3 -g -funroll-loops -fprefetch-loop-arrays -fpermissive -fno-exceptions -static-libgcc -Wl,--hash-style=both,--as-needed -DPARSEC_VERSION=3.0-beta-20150206 -DENABLE_THREADS -pthread -c streamcluster.cpp


g++ -pg -O3 -g -funroll-loops -fprefetch-loop-arrays -fpermissive -fno-exceptions -static-libgcc -Wl,--hash-style=both,--as-needed -DPARSEC_VERSION=3.0-beta-20150206 -DENABLE_THREADS -pthread -c parsec_barrier.cpp

g++ -pg -O3 -g -funroll-loops -fprefetch-loop-arrays -fpermissive -fno-exceptions -static-libgcc -Wl,--hash-style=both,--as-needed -DPARSEC_VERSION=3.0-beta-20150206 -DENABLE_THREADS -pthread -L/usr/lib64 -L/usr/lib streamcluster.o parsec_barrier.o  -o streamcluster

#### Parsec
parsecmgmt -a run -p streamcluster -i native

#### Geral
./streamcluster 10 20 128 1000000 200000 5000 none output.txt 4
./streamcluster 2 5 1 10 10 5 none output.txt 2 

gprof streamcluster gmon.out

gprof streamcluster | gprof2dot > out.dot

sudo apt-get install -y xdot

xdot out.dot


#### Pascal

pascalanalyzer ./streamcluster -c 1:4 -i input.txt -o main.json

pascalanalyzer -c 1:4 -i "10 20 128 1000000 200000 5000 none output.txt __nt__" ./streamcluster -o main.json

#### Open mp

g++ -pg -O3 -g -funroll-loops -Wall -fopenmp -fprefetch-loop-arrays -fpermissive -fno-exceptions -static-libgcc -Wl,--hash-style=both,--as-needed -DPARSEC_VERSION=3.0-beta-20150206 -c streamcluster.cpp
g++ -pg -O3 -g -funroll-loops -Wall -fopenmp -fprefetch-loop-arrays -fpermissive -fno-exceptions -static-libgcc -Wl,--hash-style=both,--as-needed -DPARSEC_VERSION=3.0-beta-20150206 -L/usr/lib64 -L/usr/lib streamcluster.o -o streamcluster

#### Serial

g++ -pg -O3 -g -funroll-loops -fprefetch-loop-arrays -fpermissive -fno-exceptions -static-libgcc -Wl,--hash-style=both,--as-needed -DPARSEC_VERSION=3.0-beta-20150206 -c streamcluster.cpp
g++ -pg -O3 -g -funroll-loops -fprefetch-loop-arrays -fpermissive -fno-exceptions -static-libgcc -Wl,--hash-style=both,--as-needed -DPARSEC_VERSION=3.0-beta-20150206 -L/usr/lib64 -L/usr/lib streamcluster.o -o streamcluster

#### Pthreads - OMP
g++ -pg -O3 -g -funroll-loops -Wall -fopenmp -fprefetch-loop-arrays -fpermissive -fno-exceptions -static-libgcc -Wl,--hash-style=both,--as-needed -DPARSEC_VERSION=3.0-beta-20150206 -DENABLE_THREADS -pthread -c streamcluster.cpp

g++ -pg -O3 -g -funroll-loops -Wall -fopenmp -fprefetch-loop-arrays -fpermissive -fno-exceptions -static-libgcc -Wl,--hash-style=both,--as-needed -DPARSEC_VERSION=3.0-beta-20150206 -DENABLE_THREADS -pthread -c parsec_barrier.cpp

g++ -pg -O3 -g -funroll-loops -Wall -fopenmp -fprefetch-loop-arrays -fpermissive -fno-exceptions -static-libgcc -Wl,--hash-style=both,--as-needed -DPARSEC_VERSION=3.0-beta-20150206 -DENABLE_THREADS -pthread -L/usr/lib64 -L/usr/lib streamcluster.o parsec_barrier.o  -o streamcluster
