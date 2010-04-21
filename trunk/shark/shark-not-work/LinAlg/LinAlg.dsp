# Microsoft Developer Studio Project File - Name="LinAlg" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** NICHT BEARBEITEN **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=LinAlg - Win32 Debug
!MESSAGE Dies ist kein g�ltiges Makefile. Zum Erstellen dieses Projekts mit NMAKE
!MESSAGE verwenden Sie den Befehl "Makefile exportieren" und f�hren Sie den Befehl
!MESSAGE 
!MESSAGE NMAKE /f "LinAlg.mak".
!MESSAGE 
!MESSAGE Sie k�nnen beim Ausf�hren von NMAKE eine Konfiguration angeben
!MESSAGE durch Definieren des Makros CFG in der Befehlszeile. Zum Beispiel:
!MESSAGE 
!MESSAGE NMAKE /f "LinAlg.mak" CFG="LinAlg - Win32 Debug"
!MESSAGE 
!MESSAGE F�r die Konfiguration stehen zur Auswahl:
!MESSAGE 
!MESSAGE "LinAlg - Win32 Release" (basierend auf  "Win32 (x86) Static Library")
!MESSAGE "LinAlg - Win32 Debug" (basierend auf  "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "LinAlg - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\lib\WinNT\Release"
# PROP Intermediate_Dir "..\win_tmp\Release\LinAlg"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /w /W0 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /w /W0 /GR /GX /O2 /I "../include" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "LinAlg - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\lib\WinNT\Debug"
# PROP Intermediate_Dir "..\win_tmp\Debug\LinAlg"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /w /W0 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MD /w /W0 /GR /GX /Z7 /Od /I "../include" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "LinAlg - Win32 Release"
# Name "LinAlg - Win32 Debug"
# Begin Group "Quellcodedateien"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\src\bfgs.cpp
# End Source File
# Begin Source File

SOURCE=.\src\bfgs2.cpp
# End Source File
# Begin Source File

SOURCE=.\src\cblnsrch.cpp
# End Source File
# Begin Source File

SOURCE=.\src\detsymm.cpp
# End Source File
# Begin Source File

SOURCE=.\src\discrimAnalysis.cpp
# End Source File
# Begin Source File

SOURCE=.\src\dlinmin.cpp
# End Source File
# Begin Source File

SOURCE=.\src\eigenerr.cpp
# End Source File
# Begin Source File

SOURCE=.\src\eigensort.cpp
# End Source File
# Begin Source File

SOURCE=.\src\eigensymm.cpp
# End Source File
# Begin Source File

SOURCE=.\src\eigensymmJacobi.cpp
# End Source File
# Begin Source File

SOURCE=.\src\eigensymmJacobi2.cpp
# End Source File
# Begin Source File

SOURCE=.\src\fft.cpp
# End Source File
# Begin Source File

SOURCE=.\src\g_inverse.cpp
# End Source File
# Begin Source File

SOURCE=.\src\invert.cpp
# End Source File
# Begin Source File

SOURCE=.\src\linalg.cpp
# End Source File
# Begin Source File

SOURCE=.\src\LinearClassifier.cpp
# End Source File
# Begin Source File

SOURCE=.\src\linearRegress.cpp
# End Source File
# Begin Source File

SOURCE=.\src\LinearRegression.cpp
# End Source File
# Begin Source File

SOURCE=.\src\linmin.cpp
# End Source File
# Begin Source File

SOURCE=.\src\lnsrch.cpp
# End Source File
# Begin Source File

SOURCE=.\src\PCA.cpp
# End Source File
# Begin Source File

SOURCE=.\src\rank.cpp
# End Source File
# Begin Source File

SOURCE=.\src\rankDecomp.cpp
# End Source File
# Begin Source File

SOURCE=.\src\svd.cpp
# End Source File
# Begin Source File

SOURCE=.\src\svdrank.cpp
# End Source File
# Begin Source File

SOURCE=.\src\svdsort.cpp
# End Source File
# End Group
# Begin Group "Header-Dateien"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\include\LinAlg\arrayoptimize.h
# End Source File
# Begin Source File

SOURCE=..\include\LinAlg\fft.h
# End Source File
# Begin Source File

SOURCE=..\include\LinAlg\linalg.h
# End Source File
# Begin Source File

SOURCE=..\include\LinAlg\LinearClassifier.h
# End Source File
# Begin Source File

SOURCE=..\include\LinAlg\LinearRegression.h
# End Source File
# Begin Source File

SOURCE=..\include\LinAlg\PCA.h
# End Source File
# End Group
# End Target
# End Project