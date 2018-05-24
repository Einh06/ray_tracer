@echo off

SET OPTS=-FC -GR- -EHa- -nologo -Zi /O2
SET CODE_HOME=%cd%
pushd build\
cl %OPTS% %CODE_HOME%\main.c -Femain
popd
