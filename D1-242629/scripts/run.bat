@echo off
REM run.bat - Windows batch file to run spMV_parallel
REM This script is in: project/scripts/

REM Get the script's directory and go to project root
cd /d "%~dp0.."

REM Check if matrix file is provided
if "%~1"=="" (
    echo Usage: scripts\run.bat ^<matrix.csr^>
    echo Example: scripts\run.bat matrices\m1.csr
    exit /b 1
)

set MATRIX=%~1
set EXE=src\spMV_parall.exe

REM Check if executable exists
if not exist "%EXE%" (
    echo Error: Executable not found at %EXE%
    echo Please compile first with: gcc -fopenmp -o src/spMV_parall src/spMV_parall.c
    exit /b 1
)

REM Check if matrix file exists
if not exist "%MATRIX%" (
    echo Error: Matrix file not found: %MATRIX%
    exit /b 1
)

echo ===================================
echo Running SpMV OpenMP Benchmarks
echo Matrix: %MATRIX%
echo ===================================
echo.

REM Test different thread counts with static schedule
echo --- Testing thread scaling (static, chunk=100) ---
%EXE% 1 static 100 "%MATRIX%"
%EXE% 2 static 100 "%MATRIX%"
%EXE% 4 static 100 "%MATRIX%"
%EXE% 8 static 100 "%MATRIX%"
%EXE% 16 static 100 "%MATRIX%"
%EXE% 32 static 100 "%MATRIX%"
echo.

REM Test different schedules with 4 threads
echo --- Testing schedules (4 threads, chunk=100) ---
%EXE% 4 static 100 "%MATRIX%"
%EXE% 4 dynamic 100 "%MATRIX%"
%EXE% 4 guided 100 "%MATRIX%"
%EXE% 4 auto 0 "%MATRIX%"
echo.

REM Test different chunk sizes with dynamic schedule
echo --- Testing chunk sizes (4 threads, dynamic) ---
%EXE% 4 dynamic 1 "%MATRIX%"
%EXE% 4 dynamic 10 "%MATRIX%"
%EXE% 4 dynamic 100 "%MATRIX%"
%EXE% 4 dynamic 1000 "%MATRIX%"
%EXE% 4 dynamic 0 "%MATRIX%"
echo.

echo ===================================
echo Benchmarks Complete
echo ===================================