#!/bin/bash

cd src

make test

./test | python visualize.py

cd ..