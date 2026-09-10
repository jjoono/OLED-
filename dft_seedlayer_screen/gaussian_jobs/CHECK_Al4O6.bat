@echo off
rem  Al4O6 폴더의 계산 결과가 물리적으로 말이 되는지 점검합니다.
rem  더블클릭하면 창에 결과가 나옵니다. 그 내용을 그대로 복사해서 보내주세요.
setlocal enabledelayedexpansion
cd /d "%~dp0Al4O6"

echo ================================================================
echo  Al4O6 결과 점검
echo ================================================================
echo.
echo  파일               정상종료  SCF 에너지 (Hartree)      스핀 S**2
echo  ----------------------------------------------------------------

for %%F in (*.out) do (
  set "NORM=NO "
  findstr /c:"Normal termination" "%%F" >nul 2>&1
  if not errorlevel 1 set "NORM=YES"

  set "E=-"
  for /f "tokens=5" %%A in ('findstr /c:"SCF Done" "%%F" 2^>nul') do set "E=%%A"

  set "S2=-"
  for /f "tokens=5" %%B in ('findstr /c:"S**2 before annihilation" "%%F" 2^>nul') do set "S2=%%B"

  set "NAME=%%~nF                    "
  echo  !NAME:~0,18! !NORM!      !E!        !S2!
)

echo.
echo  ----------------------------------------------------------------
echo  기대값:
echo    정상종료   : 8개 전부 YES
echo    SCF 에너지 : 전부 -1565 근처, 서로 0.02 이내
echo    S**2       : 0.75 근처 (홑전자 1개인 이중항)
echo.
echo  위 표를 통째로 복사해서 보내주시면 장벽을 계산해 드립니다.
echo ================================================================
pause
