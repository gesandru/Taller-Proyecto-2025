# Taller-Proyecto-2025

Para compilar los tests de matrix:
g++ tests/tests_matrix.cpp src/*.cpp -lm -std=c++23 -o bin/tests_matrix.exe

Para compilar los tests unitarios del resto de funciones:
g++ tests/tests.cpp src/*.cpp -lm -std=c++23 -o bin/tests.exe


Para compilar el test de la aplicación completa:
g++ tests/EKF_GEOS3.cpp src/*cpp -lm -std=c++23 -o bin/main.exe

Documentación en carpeta doc, acceder al index.html