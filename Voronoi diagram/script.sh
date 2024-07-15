#!/bin/bash

g++-13 main.cpp
./a.out < input.txt > results.txt
python3 plot.py