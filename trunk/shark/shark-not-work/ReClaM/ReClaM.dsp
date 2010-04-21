# Microsoft Developer Studio Project File - Name="ReClaM" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=ReClaM - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "ReClaM.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "ReClaM.mak" CFG="ReClaM - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "ReClaM - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "ReClaM - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "ReClaM - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\lib\WinNT\Release"
# PROP Intermediate_Dir "..\win_tmp\Release\ReClaM"
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

!ELSEIF  "$(CFG)" == "ReClaM - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\lib\WinNT\Debug"
# PROP Intermediate_Dir "..\win_tmp\Debug\ReClaM"
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

# Name "ReClaM - Win32 Release"
# Name "ReClaM - Win32 Debug"
# Begin Group "Quellcodedateien"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\src\BFGS.cpp
# End Source File
# Begin Source File

SOURCE=.\src\CG.cpp
# End Source File
# Begin Source File

SOURCE=.\src\EarlyStopping.cpp
# End Source File
# Begin Source File

SOURCE=.\src\FFNet.cpp
# End Source File
# Begin Source File

SOURCE=.\src\FileUtil.cpp
# End Source File
# Begin Source File

SOURCE=.\src\MSEFFNet.cpp
# End Source File
# Begin Source File

SOURCE=.\src\MSERNNet.cpp
# End Source File
# Begin Source File

SOURCE=.\src\NetParams.cpp
# End Source File
# Begin Source File

SOURCE=.\src\Params.cpp
# End Source File
# Begin Source File

SOURCE=.\src\VarianceEstimator.cpp
# End Source File
# End Group
# Begin Group "Header-Dateien"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\include\ReClaM\AdpBP.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\BFGS.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\CE.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\CG.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\CMAC.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\createConnectionMatrix.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\CrossEntropy.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\CrossEntropyIndependent.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\DF_CrossEntropy.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\DF_MeanSquaredError.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\EarlyStopping.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\EnhancedRprop.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\ErrorMeasures.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\FFNet.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\FFNets.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\FFNetSource.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\FileUtil.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\IOTools.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\LinOutFFNet.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\LMSEFFNet.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\MeanSquaredError.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\ModelInterface.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\MSEBFFNet.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\MSEFFNet.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\MSERBFNet.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\MSERNNet.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\NetParams.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\Paraboloid.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\Params.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\ProbNet.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\Quickprop.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\RBFNet.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\Rprop.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\SquaredError.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\SteepestDescent.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\StochasticGradientDescent.h
# End Source File
# Begin Source File

SOURCE=..\include\ReClaM\VarianceEstimator.h
# End Source File
# End Group
# End Target
# End Project