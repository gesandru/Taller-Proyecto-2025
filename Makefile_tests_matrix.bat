g++ tests/tests_matrix.cpp src/*.cpp -lm -std=c++23 -o bin/tests_matrix.exe
cd bin
tests_matrix.exe
cd ..
pause