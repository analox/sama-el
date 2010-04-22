# Microsoft Developer Studio Project File - Name="TimeSeries" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** NICHT BEARBEITEN **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=TimeSeries - Win32 Debug
!MESSAGE Dies ist kein g�ltiges Makefile. Zum Erstellen dieses Projekts mit NMAKE
!MESSAGE verwenden Sie den Befehl "Makefile exportieren" und f�hren Sie den Befehl
!MESSAGE 
!MESSAGE NMAKE /f "TimeSeries.mak".
!MESSAGE 
!MESSAGE Sie k�nnen beim Ausf�hren von NMAKE eine Konfiguration angeben
!MESSAGE durch Definieren des Makros CFG in der Befehlszeile. Zum Beispiel:
!MESSAGE 
!MESSAGE NMAKE /f "TimeSeries.mak" CFG="TimeSeries - Win32 Debug"
!MESSAGE 
!MESSAGE F�r die Konfiguration stehen zur Auswahl:
!MESSAGE 
!MESSAGE "TimeSeries - Win32 Release" (basierend auf  "Win32 (x86) Static Library")
!MESSAGE "TimeSeries - Win32 Debug" (basierend auf  "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "TimeSeries - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\..\lib\WinNT\Release"
# PROP Intermediate_Dir "..\..\win_tmp\Release\TimeSeries"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GR /GX /O2 /I "..\..\include" /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD BASE RSC /l 0x407 /d "NDEBUG"
# ADD RSC /l 0x407 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "TimeSeries - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\..\lib\WinNT\Debug"
# PROP Intermediate_Dir "..\..\win_tmp\Debug\TimeSeries"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MD /W3 /GR /GX /Z7 /Od /I "..\..\include" /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
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

# Name "TimeSeries - Win32 Release"
# Name "TimeSeries - Win32 Debug"
# Begin Group "Quellcodedateien"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\src\BimodalBrownianProcess.cpp
# End Source File
# Begin Source File

SOURCE=.\src\DiscreteMackeyGlass.cpp
# End Source File
# Begin Source File

SOURCE=.\src\Lorenz63.cpp
# End Source File
# Begin Source File

SOURCE=.\src\Lorenz84.cpp
# End Source File
# Begin Source File

SOURCE=.\src\MackeyGlass.cpp
# End Source File
# Begin Source File

SOURCE=.\src\NoisyIOSamples.cpp
# End Source File
# Begin Source File

SOURCE=.\src\NoisyMackeyGlass.cpp
# End Source File
# Begin Source File

SOURCE=".\src\RK4-1D.cpp"
# End Source File
# Begin Source File

SOURCE=.\src\RK4.cpp
# End Source File
# End Group
# Begin Group "Header-Dateien"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\..\include\TestData\TimeSeries\BimodalBrownianProcess.h
# End Source File
# Begin Source File

SOURCE=..\..\include\TestData\TimeSeries\Counter.h
# End Source File
# Begin Source File

SOURCE=..\..\include\TestData\TimeSeries\DiscreteMackeyGlass.h
# End Source File
# Begin Source File

SOURCE=..\..\include\TestData\TimeSeries\Embedding.h
# End Source File
# Begin Source File

SOURCE=..\..\include\TestData\TimeSeries\Generator.h
# End Source File
# Begin Source File

SOURCE=..\..\include\TestData\TimeSeries\IOGenerator.h
# End Source File
# Begin Source File

SOURCE=..\..\include\TestData\TimeSeries\IOSamples.h
# End Source File
# Begin Source File

SOURCE=..\..\include\TestData\TimeSeries\Lorenz63.h
# End Source File
# Begin Source File

SOURCE=..\..\include\TestData\TimeSeries\Lorenz84.h
# End Source File
# Begin Source File

SOURCE=..\..\include\TestData\TimeSeries\MackeyGlass.h
# End Source File
# Begin Source File

SOURCE=..\..\include\TestData\TimeSeries\NoisyIOSamples.h
# End Source File
# Begin Source File

SOURCE=..\..\include\TestData\TimeSeries\NoisyMackeyGlass.h
# End Source File
# Begin Source File

SOURCE="..\..\include\TestData\TimeSeries\RK4-1D.h"
# End Source File
# Begin Source File

SOURCE=..\..\include\TestData\TimeSeries\RK4.h
# End Source File
# Begin Source File

SOURCE=..\..\include\TestData\TimeSeries\SelectComponent.h
# End Source File
# End Group
# End Target
# End Project