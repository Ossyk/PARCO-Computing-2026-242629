@echo off
:: main_menu_windows.bat - Interactive menu for SpMV project
:: Author: 242629
:: This script should be in project\scripts\

cd /d "%~dp0\.."

echo Welcome to 242629's script!
echo.

:menu
echo ==========================================
echo === Main Menu ===
echo 1^) Create your matrix
echo 2^) Run sequential code
echo 3^) Run parallel code
echo 4^) Run default benchmark code
echo 5^) Run personalized benchmark code
echo 6^) See available matrices
echo 7^) Exit
echo ==========================================
echo.

set /p choice="Enter your choice (1-7): "
echo.

if "%choice%"=="1" goto create_matrix
if "%choice%"=="2" goto run_seq
if "%choice%"=="3" goto run_par
if "%choice%"=="4" goto run_default
if "%choice%"=="5" goto run_personalized
if "%choice%"=="6" goto list_matrices
if "%choice%"=="7" goto goodbye

echo Invalid choice. Please try again.
echo.
goto menu


:create_matrix
echo Creating your square CSR matrix...
set /p size=Enter the number of rows and columns:
set /p sparsity=Enter sparsity:
set /p name=Enter matrix name (other than m1-m5 and already chosen names):

gcc src\matrix_generator.c src\csr_utils.c -o src\matrix_generator.exe
if exist src\matrix_generator.exe (
    src\matrix_generator.exe %size% %size% %sparsity% matrices\%name%.csr
)
echo.
goto menu


:run_seq
echo Running sequential code...
set /p name=Enter matrix name (from m1 to m5 or created):

gcc src\spMV_seq.c src\csr_utils.c -o src\spMV_seq.exe
if exist src\spMV_seq.exe (
    src\spMV_seq.exe matrices\%name%.csr
)
echo.
goto menu


:run_par
echo Running parallel code...
set /p threads=Enter number of threads:
set /p scheduler=Enter scheduler:
set /p chunks=Enter chunk size:
set /p name=Enter matrix name:

gcc -fopenmp src\spMV_parall.c src\csr_utils.c -o src\spMV_parall.exe
if exist src\spMV_parall.exe (
    src\spMV_parall.exe %threads% %scheduler% %chunks% matrices\%name%.csr
)
echo.
goto menu


:run_default
echo Running default benchmark test...
set /p name=Enter matrix name (from m1 to m5 or created):
if exist scripts\run_windows.bat (
    call scripts\run_windows.bat 4 dynamic 100 matrices\%name%.csr
)
echo.
goto menu


:run_personalized
echo Running personalized benchmark test...
set /p threads=Enter number of threads:
set /p scheduler=Enter scheduler:
set /p chunks=Enter chunk size:
set /p name=Enter matrix name:

if exist scripts\run_windows.bat (
    call scripts\run_windows.bat %threads% %scheduler% %chunks% matrices\%name%.csr
)
echo.
goto menu


:list_matrices
echo Available matrices:
dir /b matrices
echo.
goto menu


:goodbye
echo Goodbye!
exit
