@echo off
rem ---------------------------------------------------------------------------
rem Run the E_d campaign with Gaussian 16W, without depending on this machine's
rem PATH being intact.
rem
rem Gaussian is not one program: g16.exe launches a chain of link executables
rem (l101.exe, l502.exe, ...) that it finds through GAUSS_EXEDIR and PATH. On a
rem machine whose PATH has been damaged, g16 starts and dies immediately with no
rem useful message, so both are set here explicitly rather than assumed.
rem
rem   run_gaussian_windows.bat --selftest    check input and parser, no Gaussian
rem   run_gaussian_windows.bat pyridine      one system
rem   run_gaussian_windows.bat               all 26
rem
rem Edit the four SET lines below if anything sits somewhere else.
rem ---------------------------------------------------------------------------
setlocal

set "G16DIR=C:\G16W"
set "PY=C:\Users\JHKIM\miniforge3\python.exe"
set "GAUSS_NPROC=8"
set "GAUSS_MEM_GB=16"

rem --- rebuild a working PATH from scratch, so a damaged one cannot matter ----
set "PATH=%SystemRoot%\System32;%SystemRoot%;%SystemRoot%\System32\Wbem;%G16DIR%;%PATH%"
set "GAUSS_EXEDIR=%G16DIR%"
set "GAUSS_EXE=%G16DIR%\g16.exe"
if not defined GAUSS_SCRDIR set "GAUSS_SCRDIR=%TEMP%"

rem --- check before running, so a missing piece is named ----------------------
if not exist "%GAUSS_EXE%" (
  echo.
  echo   Gaussian not found at %GAUSS_EXE%
  echo   Edit G16DIR at the top of this file.
  echo.
  exit /b 1
)
if not exist "%PY%" (
  echo.
  echo   Python not found at %PY%
  echo.
  echo   Install Miniforge ^(it brings numpy with it^):
  echo     https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe
  echo   Choose "Just Me" and keep the default location, then run this again.
  echo.
  echo   If Python is already somewhere else, edit PY at the top of this file.
  echo.
  exit /b 1
)
"%PY%" -c "import numpy" 2>nul
if errorlevel 1 (
  echo installing numpy...
  "%PY%" -m pip install numpy
  if errorlevel 1 (
    echo.
    echo   numpy could not be installed. If this machine has no internet,
    echo   copy a numpy wheel over and: "%PY%" -m pip install numpy-....whl
    echo.
    exit /b 1
  )
)

cd /d "%~dp0.."
echo python  : %PY%
echo gaussian: %GAUSS_EXE%
echo cores   : %GAUSS_NPROC%    memory: %GAUSS_MEM_GB% GB
echo scratch : %GAUSS_SCRDIR%
echo.

"%PY%" scripts\102_ed_gaussian.py %*
