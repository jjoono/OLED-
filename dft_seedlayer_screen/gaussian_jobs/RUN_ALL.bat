@echo off
rem  Start the campaign. Everything real happens in run_campaign.py -- this file
rem  only finds a Python and hands over, because nested batch loops cannot run
rem  several Gaussian jobs at once or keep a job order that is not alphabetical.
cd /d "%~dp0"

where python >nul 2>&1
if not errorlevel 1 (
  python run_campaign.py %*
  goto :done
)
if exist "%USERPROFILE%\miniforge3\python.exe" (
  "%USERPROFILE%\miniforge3\python.exe" run_campaign.py %*
  goto :done
)
echo   No Python found on PATH.
echo   Open the Miniforge Prompt and run:  python run_campaign.py

:done
echo.
pause
