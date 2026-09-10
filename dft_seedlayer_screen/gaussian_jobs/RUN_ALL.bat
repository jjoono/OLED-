@echo off
rem  Run every Gaussian job in this folder tree, in order.
rem  Safe to stop and restart: a job whose .log already exists is skipped.
setlocal enabledelayedexpansion

set "G16DIR=C:\G16W"
set "PATH=%SystemRoot%\System32;%SystemRoot%;%G16DIR%;%PATH%"
set "GAUSS_EXEDIR=%G16DIR%"
if not defined GAUSS_SCRDIR set "GAUSS_SCRDIR=%TEMP%"

if not exist "%G16DIR%\g16.exe" (
  echo   Gaussian not found at %G16DIR%\g16.exe
  echo   Edit G16DIR at the top of this file.
  pause & exit /b 1
)

cd /d "%~dp0"
set /a done=0
set /a skip=0

for /r %%F in (*.gjf) do (
  if exist "%%~dpnF.log" (
    set /a skip+=1
  ) else (
    echo [!time!] %%~nF
    pushd "%%~dpF"
    "%G16DIR%\g16.exe" "%%~nxF"
    popd
    set /a done+=1
  )
)

echo.
echo   finished: !done! run, !skip! already had logs
echo   Now send the whole folder back, or just the .log files.
pause
