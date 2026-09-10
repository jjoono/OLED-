@echo off
rem  Summarise every folder that has finished. Safe to run at any time.
rem  ASCII only on purpose: a .bat with non-ASCII text is read by cmd in the
rem  console code page, not UTF-8, and comes out as mojibake.
setlocal enabledelayedexpansion
cd /d "%~dp0"

echo ================================================================
echo  Gaussian results so far
echo ================================================================
echo.

set /a NDONE=0
set /a NTODO=0

for /d %%D in (*) do (
  set /a N=0
  set /a OK=0
  for %%F in ("%%D\*.gjf") do set /a N+=1
  for %%F in ("%%D\*.out") do (
    findstr /c:"Normal termination" "%%F" >nul 2>&1
    if not errorlevel 1 set /a OK+=1
  )
  if !N! gtr 0 (
    if !OK! equ !N! (
      set /a NDONE+=1
      echo [DONE] %%D   !OK!/!N!
      for %%F in ("%%D\*.out") do (
        set "E=-"
        for /f "tokens=5" %%A in ('findstr /c:"SCF Done" "%%F" 2^>nul') do set "E=%%A"
        set "S2=-"
        for /f "tokens=4" %%B in ('findstr /c:"S**2 before annihilation" "%%F" 2^>nul') do set "S2=%%B"
        set "NAME=%%~nF                        "
        echo        !NAME:~0,22! !E!   S2=!S2!
      )
      echo.
    ) else (
      set /a NTODO+=1
      echo [....] %%D   !OK!/!N!
    )
  )
)

echo ================================================================
echo   !NDONE! folders complete, !NTODO! still running or waiting
echo.
echo   Copy everything above and send it back.
echo   S2 should be near 0.75. Far from it means the SCF found a
echo   spin-contaminated state and that energy is not the doublet's.
echo ================================================================
pause
