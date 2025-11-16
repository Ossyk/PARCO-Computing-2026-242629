@echo off

cd /d "%~dp0.."

set THREADS=%1
set SCHEDULER=%2
set CHUNKS=%3
set MATRIX=%4
set EXE=src\spMV_parall.exe

:: Check if executable exists
if not exist "%EXE%" (
    echo Error: Executable not found at %EXE%
    echo Please compile first with:
    echo    gcc -fopenmp -o src\spMV_parall.exe src\spMV_parall.c
    pause
    exit /b 1
)

:: Check if matrix file exists
if not exist "%MATRIX%" (
    echo Error: Matrix file not found: %MATRIX%
    pause
    exit /b 1
)

echo ===================================
echo Running SpMV OpenMP Benchmarks
echo Matrix: %MATRIX%
echo ===================================
echo.

:: Test different thread scaling
echo --- Testing thread scaling ---
%EXE% 1 %SCHEDULER% %CHUNKS% %MATRIX%
%EXE% 2 %SCHEDULER% %CHUNKS% %MATRIX%
%EXE% 4 %SCHEDULER% %CHUNKS% %MATRIX%
%EXE% 8 %SCHEDULER% %CHUNKS% %MATRIX%
%EXE% 16 %SCHEDULER% %CHUNKS% %MATRIX%
%EXE% 32 %SCHEDULER% %CHUNKS% %MATRIX%
%EXE% 64 %SCHEDULER% %CHUNKS% %MATRIX%

echo.

:: Test different schedules 
echo --- Testing schedules ---
%EXE% %THREADS% static %CHUNKS% %MATRIX%
%EXE% %THREADS% dynamic %CHUNKS% %MATRIX%
%EXE% %THREADS% guided %CHUNKS% %MATRIX%
%EXE% %THREADS% auto 0 %MATRIX%
echo.

:: Test different chunk sizes
echo --- Testing chunk sizes ---
%EXE% %THREADS% %SCHEDULER% 1 %MATRIX%
%EXE% %THREADS% %SCHEDULER% 10 %MATRIX%
%EXE% %THREADS% %SCHEDULER% 100 %MATRIX%
%EXE% %THREADS% %SCHEDULER% 1000 %MATRIX%
%EXE% 4 dynamic 0 %MATRIX%
echo.

echo ===================================
echo Benchmarks Complete
echo ===================================

pause


