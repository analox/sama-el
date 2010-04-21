# Microsoft Developer Studio Project File - Name="MOOEALib" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=MOOEALib - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "MOOEALib.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "MOOEALib.mak" CFG="MOOEALib - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "MOOEALib - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "MOOEALib - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "MOOEALib - Win32 Release"

# PROP BASE Use_MFC 2
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\lib\WinNT\Release"
# PROP Intermediate_Dir "..\win_tmp\Release\MOO-EALib"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MD /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_AFXDLL" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /I "..\include" /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "__NO_GENERIC_IOSTREAM" /Fp"..\win_tmp\Release\MOO-EALib/MOO-EALib.pch" /YX /FD /c
# ADD BASE RSC /l 0x407 /d "NDEBUG" /d "_AFXDLL"
# ADD RSC /l 0x407 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo /o"..\lib\WinNT\Release/MOO-EALib.bsc"
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\lib\WinNT\Release\MOO-EALib.lib"

!ELSEIF  "$(CFG)" == "MOOEALib - Win32 Debug"

# PROP BASE Use_MFC 2
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\lib\WinNT\Debug"
# PROP Intermediate_Dir "..\win_tmp\Debug\MOO-EALib"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /MDd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_AFXDLL" /YX /FD /GZ /c
# ADD CPP /nologo /MD /W3 /GR /GX /Z7 /Od /I "../include" /D "_DEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "__NO_GENERIC_IOSTREAM" /YX /FD /GZ /c
# ADD BASE RSC /l 0x407 /d "_DEBUG" /d "_AFXDLL"
# ADD RSC /l 0x407 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo /o"..\lib\WinNT\Debug/MOO-EALib.bsc"
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\lib\WinNT\Debug\MOO-EALib.lib"

!ENDIF 

# Begin Target

# Name "MOOEALib - Win32 Release"
# Name "MOOEALib - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\src\ArchiveMOO.cpp
# End Source File
# Begin Source File

SOURCE=.\src\IndividualMOO.cpp
# End Source File
# Begin Source File

SOURCE=.\src\PopulationMOO.cpp
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE="..\include\MOO-EALib\ArchiveMOO.h"
# End Source File
# Begin Source File

SOURCE="..\include\MOO-EALib\IndividualMOO.h"
# End Source File
# Begin Source File

SOURCE="..\include\MOO-EALib\PopulationMOO.h"
# End Source File
# Begin Source File

SOURCE="..\include\MOO-EALib\TestFunction.h"
# End Source File
# End Group
# End Target
# End Project
