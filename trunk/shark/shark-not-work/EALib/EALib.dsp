# Microsoft Developer Studio Project File - Name="EALib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=EALib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "EALib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "EALib.mak" CFG="EALib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "EALib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "EALib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "EALib - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\lib\WinNT\Release"
# PROP Intermediate_Dir "..\win_tmp\Release\EALib"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /I "..\include" /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "__NO_GENERIC_IOSTREAM" /YX /FD /c
# ADD BASE RSC /l 0x407 /d "NDEBUG"
# ADD RSC /l 0x407 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "EALib - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\lib\WinNT\Debug"
# PROP Intermediate_Dir "..\win_tmp\Debug\EALib"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MD /W3 /GR /GX /Z7 /Od /I "../include" /D "_DEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "__NO_GENERIC_IOSTREAM" /YX /FD /GZ /c
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

# Name "EALib - Win32 Release"
# Name "EALib - Win32 Debug"
# Begin Group "Quellcodedateien"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\src\Chromosome.cpp
# End Source File
# Begin Source File

SOURCE=.\src\ChromosomeT_bool.cpp
# End Source File
# Begin Source File

SOURCE=.\src\ChromosomeT_char.cpp
# End Source File
# Begin Source File

SOURCE=.\src\ChromosomeT_double.cpp
# End Source File
# Begin Source File

SOURCE=.\src\ChromosomeT_int.cpp
# End Source File
# Begin Source File

SOURCE=.\src\Individual.cpp
# End Source File
# Begin Source File

SOURCE=.\src\Population.cpp
# End Source File
# Begin Source File

SOURCE=.\src\PVMinterface.cpp
# End Source File
# End Group
# Begin Group "Header-Dateien"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\include\EALib\Chromosome.h
# End Source File
# Begin Source File

SOURCE=..\include\EALib\ChromosomeT.h
# End Source File
# Begin Source File

SOURCE=..\include\EALib\cma.h
# End Source File
# Begin Source File

SOURCE=..\include\EALib\cma2004.h
# End Source File
# Begin Source File

SOURCE=..\include\EALib\Individual.h
# End Source File
# Begin Source File

SOURCE=..\include\EALib\Interval.h
# End Source File
# Begin Source File

SOURCE=..\include\EALib\Population.h
# End Source File
# Begin Source File

SOURCE=..\include\EALib\PVMinterface.h
# End Source File
# Begin Source File

SOURCE=..\include\EALib\sqr.h
# End Source File
# End Group
# End Target
# End Project
