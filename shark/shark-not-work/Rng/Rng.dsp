# Microsoft Developer Studio Project File - Name="Rng" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** NICHT BEARBEITEN **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=Rng - Win32 Debug
!MESSAGE Dies ist kein g�ltiges Makefile. Zum Erstellen dieses Projekts mit NMAKE
!MESSAGE verwenden Sie den Befehl "Makefile exportieren" und f�hren Sie den Befehl
!MESSAGE 
!MESSAGE NMAKE /f "Rng.mak".
!MESSAGE 
!MESSAGE Sie k�nnen beim Ausf�hren von NMAKE eine Konfiguration angeben
!MESSAGE durch Definieren des Makros CFG in der Befehlszeile. Zum Beispiel:
!MESSAGE 
!MESSAGE NMAKE /f "Rng.mak" CFG="Rng - Win32 Debug"
!MESSAGE 
!MESSAGE F�r die Konfiguration stehen zur Auswahl:
!MESSAGE 
!MESSAGE "Rng - Win32 Release" (basierend auf  "Win32 (x86) Static Library")
!MESSAGE "Rng - Win32 Debug" (basierend auf  "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "Rng - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\lib\WinNT\Release"
# PROP Intermediate_Dir "..\win_tmp\Release\Rng"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /I "..\include" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x407 /d "NDEBUG"
# ADD RSC /l 0x407 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "Rng - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\lib\WinNT\Debug"
# PROP Intermediate_Dir "..\win_tmp\Debug\Rng"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MD /W3 /GR /GX /Z7 /Od /I "..\include" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD BASE RSC /l 0x407 /d "_DEBUG"
# ADD RSC /l 0x407 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "Rng - Win32 Release"
# Name "Rng - Win32 Debug"
# Begin Group "Quellcodedateien"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\src\Bernoulli.cpp
# End Source File
# Begin Source File

SOURCE=.\src\Binomial.cpp
# End Source File
# Begin Source File

SOURCE=.\src\Cauchy.cpp
# End Source File
# Begin Source File

SOURCE=.\src\DiffGeometric.cpp
# End Source File
# Begin Source File

SOURCE=.\src\DiscreteUniform.cpp
# End Source File
# Begin Source File

SOURCE=.\src\Erlang.cpp
# End Source File
# Begin Source File

SOURCE=.\src\Geometric.cpp
# End Source File
# Begin Source File

SOURCE=.\src\GlobalRng.cpp
# End Source File
# Begin Source File

SOURCE=.\src\HyperGeometric.cpp
# End Source File
# Begin Source File

SOURCE=.\src\LogNormal.cpp
# End Source File
# Begin Source File

SOURCE=.\src\NegExponential.cpp
# End Source File
# Begin Source File

SOURCE=.\src\Normal.cpp
# End Source File
# Begin Source File

SOURCE=.\src\Poisson.cpp
# End Source File
# Begin Source File

SOURCE=.\src\RNG.cpp
# End Source File
# Begin Source File

SOURCE=.\src\Uniform.cpp
# End Source File
# Begin Source File

SOURCE=.\src\Weibull.cpp
# End Source File
# End Group
# Begin Group "Header-Dateien"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\include\Rng\Bernoulli.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\Binomial.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\Cauchy.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\DiffGeometric.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\DiscreteUniform.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\Erlang.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\Geometric.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\GlobalRng.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\HyperGeometric.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\LogNormal.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\NegExponential.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\Normal.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\Poisson.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\RandomVar.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\RNG.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\Uniform.h
# End Source File
# Begin Source File

SOURCE=..\include\Rng\Weibull.h
# End Source File
# End Group
# End Target
# End Project