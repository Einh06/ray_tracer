@echo off

SET OPTS=-FC -GR- -EHa- -nologo -Zi
SET CODE_HOME=%cd%
pushd build\
cl %OPTS% %CODE_HOME%\main.c -Femain
popd
