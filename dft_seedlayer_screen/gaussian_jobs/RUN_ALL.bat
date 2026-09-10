@echo off
rem  Run every Gaussian job in this folder tree, in order.
rem
rem  The output file is named explicitly rather than left to Gaussian.
rem  G16W on Windows writes .out, not .log, and a skip check written
rem  against .log never fires -- every job would be re-run on restart.
rem
rem  A job counts as done only if its output ends in Normal termination,
rem  so a crash is retried on the next pass instead of being skipped
rem  forever.
setlocal enabledelayedexpansion

set "G16DIR=C:\G16W"
set "PATH=%SystemRoot%\System32;%SystemRoot%;%G16DIR%;%PATH%"
set "GAUSS_EXEDIR=%G16DIR%"
set "GAUSS_SCRDIR=%TEMP%\gauscr"
if not exist "%GAUSS_SCRDIR%" mkdir "%GAUSS_SCRDIR%"

if not exist "%G16DIR%\g16.exe" (
  echo   Gaussian not found at %G16DIR%\g16.exe
  echo   Edit G16DIR at the top of this file.
  pause & exit /b 1
)

cd /d "%~dp0"
set /a done=0
set /a skip=0
set /a bad=0

for /r %%F in (*.gjf) do (
  set "OUT=%%~dpnF.out"
  set "SKIPME="
  if exist "!OUT!" (
    findstr /c:"Normal termination" "!OUT!" >nul 2>&1
    if not errorlevel 1 set "SKIPME=1"
  )
  if defined SKIPME (
    set /a skip+=1
  ) else (
    echo [!time!] %%~nF
    pushd "%%~dpF"
    "%G16DIR%\g16.exe" "%%~nxF" "%%~nF.out"
    popd
    findstr /c:"Normal termination" "!OUT!" >nul 2>&1
    if errorlevel 1 (
      set /a bad+=1
      echo        ... did NOT finish normally
    ) else (
      set /a done+=1
    )
  )
  del /q "%%~dpF\Gau-*.*" 2>nul
  del /q "%%~dpF\fort.7" 2>nul
)

echo.
echo   ran !done! successfully, !skip! already done, !bad! failed
echo.
if !bad! gtr 0 echo   Re-run this file to retry the failed ones.
echo   When all are done, send the whole gaussian_jobs folder back.
pause
