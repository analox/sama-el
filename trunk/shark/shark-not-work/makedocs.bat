REM Creates the manuals and references for all Windows-compatible
REM libraries of the Shark package.
REM You need some additional programs to get this batch file work,
REM please refer to e.g. Array\doc\README for details.

cd Array\man
CALL makeman
cd ..\doc
CALL makeref
cd ..\..\Rng\man
CALL makeman
cd ..\doc
CALL makeref
cd ..\..\EALib\man
CALL makeman
cd ..\..\LinAlg\doc
CALL makeref
cd ..\..\ReClaM\doc
call makeref
cd ..\..\EALib\doc
call makeref
cd ..\..\Metric\doc
call makeref
cd ..\..\Mixture\doc
call makeref
cd ..\..\SimAnn\doc
call makeref
cd ..\..\MOO-EALib\doc
call makeref
cd ..\..



