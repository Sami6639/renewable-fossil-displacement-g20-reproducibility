@echo off
setlocal
cd /d "%~dp0"

where conda >nul 2>nul
if errorlevel 1 (
  echo Conda was not found. Open Anaconda Prompt and run this file again.
  pause
  exit /b 1
)

call conda env list | findstr /R /C:"^g20-power-repro " >nul
if errorlevel 1 (
  echo Creating the g20-power-repro environment...
  call conda env create -f environment.yml
) else (
  echo Updating the g20-power-repro environment...
  call conda env update -n g20-power-repro -f environment.yml --prune
)
if errorlevel 1 goto :failed

echo Running the complete reproducibility workflow...
call conda run -n g20-power-repro python code\run_all.py
if errorlevel 1 goto :failed

echo.
echo SUCCESS: all reproducibility checks passed.
pause
exit /b 0

:failed
echo.
echo ERROR: the workflow did not complete. Review the message above.
pause
exit /b 1
