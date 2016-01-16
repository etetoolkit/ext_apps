#!/bin/bash

# ClustaOmega dependency
./configure --prefix $PREFIX
make
make install
make clean
