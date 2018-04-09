#!/bin/bash

CODE_HOME="$PWD"
OPTS=-g
cd build/ > /dev/null
g++ $OPTS $CODE_HOME/main -o build/
cd $CODE_HOME > /dev/null
