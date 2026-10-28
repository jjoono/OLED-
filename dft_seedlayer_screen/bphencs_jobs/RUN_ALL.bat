@echo off
rem  Start the campaign. Everything real happens in run_campaign.py -- this file
rem  only finds a Python and hands over, because nested batch loops cannot run
rem  several Gaussian jobs at once or keep a job order that is not alphabetical.
cd /d "%~dp0"

rem  A real interpreter is looked for BEFORE whatever "python" is on PATH.
rem  Windows ships a stub at
rem  %LOCALAPPDATA%\Microsoft\WindowsApps\python.exe that exists, is found by
rem  "where python", prints the word Python and does nothing. Testing that the
rem  candidate can actually run code is the only check that tells the two apart.
setlocal enabledelayedexpansion
set "PY="
for %%P in (
  "%USERPROFILE%\miniforge3\python.exe"
  "%LOCALAPPDATA%\miniforge3\python.exe"
  "%ProgramData%\miniforge3\python.exe"
  "C:\Miniforge3\python.exe"
  "%USERPROFILE%\anaconda3\python.exe"
  "%USERPROFILE%\miniconda3\python.exe"
) do if not defined PY if exist %%P set "PY=%%~P"

if not defined PY (
  for /f "delims=" %%P in ('where python 2^>nul') do (
    if not defined PY (
      echo %%P | find /i "WindowsApps" >nul || set "PY=%%P"
    )
  )
)

if defined PY (
  "!PY!" -c "import sys" >nul 2>&1
  if errorlevel 1 (
    echo   Found "!PY!" but it cannot run Python code.
    set "PY="
  )
)

if defined PY (
  echo   using !PY!
  echo.
  "!PY!" run_campaign.py %*
  goto :done
)
echo   No working Python found.
echo   The "python" on PATH is the Microsoft Store stub, which only prints the
echo   word Python. Open the Miniforge Prompt and run there:
echo       cd /d "%~dp0"
echo       python run_campaign.py

:done
echo.
pause
