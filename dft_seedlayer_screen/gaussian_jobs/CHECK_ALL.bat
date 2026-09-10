@echo off
rem  지금까지 끝난 폴더 전부의 결과를 한 번에 훑습니다.
rem  아직 안 끝난 폴더는 건너뜁니다. 언제든 돌려도 됩니다.
setlocal enabledelayedexpansion
cd /d "%~dp0"

echo ================================================================
echo  Gaussian 결과 요약
echo ================================================================
echo.

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
      echo [완료] %%D   !OK!/!N!
      for %%F in ("%%D\*.out") do (
        set "E=-"
        for /f "tokens=5" %%A in ('findstr /c:"SCF Done" "%%F" 2^>nul') do set "E=%%A"
        set "S2=-"
        for /f "tokens=5" %%B in ('findstr /c:"S**2 before annihilation" "%%F" 2^>nul') do set "S2=%%B"
        set "NAME=%%~nF                        "
        echo        !NAME:~0,22! !E!   S2=!S2!
      )
      echo.
    ) else (
      echo [진행] %%D   !OK!/!N!
    )
  )
)

echo ================================================================
echo  위 내용을 통째로 복사해서 보내주세요.
echo  S2 는 0.75 근처여야 정상입니다.
echo ================================================================
pause
