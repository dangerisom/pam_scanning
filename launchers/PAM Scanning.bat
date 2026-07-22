@echo off
REM PAM Scanning - double-click launcher (Windows).
REM
REM First run: creates the 'pam_scanning' conda environment (Python + the app +
REM NCBI BLAST+) from environment.yml, which can take a few minutes.
REM Every run after that: just opens the app.
REM
REM Requirement: install Miniforge once (a normal clickable installer):
REM   https://conda-forge.org/download/
setlocal
set ENV_NAME=pam_scanning
set HERE=%~dp0
set REPO=%HERE%..

REM Locate conda: prefer one already on PATH, else check the usual install spots.
set CONDA=
where conda >nul 2>nul && set CONDA=conda
if "%CONDA%"=="" (
  for %%C in ("%USERPROFILE%\miniforge3" "%USERPROFILE%\mambaforge" "%USERPROFILE%\miniconda3" "%USERPROFILE%\anaconda3") do (
    if exist "%%~C\Scripts\conda.exe" set CONDA="%%~C\Scripts\conda.exe"
  )
)
if "%CONDA%"=="" (
  echo Could not find conda on this PC.
  echo Please install Miniforge first ^(clickable installer^):
  echo     https://conda-forge.org/download/
  pause
  exit /b 1
)

REM Create the environment the first time only.
%CONDA% env list | findstr /b /c:"%ENV_NAME% " >nul
if errorlevel 1 (
  echo First-time setup: creating the %ENV_NAME% environment.
  echo This downloads NCBI BLAST+ and the app, and may take a few minutes...
  echo.
  pushd "%REPO%"
  %CONDA% env create -f environment.yml
  set CREATE_ERR=%errorlevel%
  popd
  if not "%CREATE_ERR%"=="0" (
    echo.
    echo Setup did not finish. Please see the messages above.
    pause
    exit /b 1
  )
  echo.
  echo Setup complete.
)

echo Starting PAM Scanning...
REM Isolate from the user's Python environment before launching.
REM
REM   * PYTHONPATH is cleared: entries there take precedence over the environment's
REM     site-packages, so any folder named pam_scanning on PYTHONPATH (e.g. an older
REM     copy of this project) shadows the installed package and the app dies with
REM     "ModuleNotFoundError: No module named 'pam_scanning.gui'".
REM   * We run from the user profile: conda run also places the current directory on
REM     sys.path, which can shadow the package the same way.
REM
REM The environment created above is self-contained, so nothing here needs PYTHONPATH.
set PYTHONPATH=
cd /d "%USERPROFILE%"
%CONDA% run --no-capture-output -n %ENV_NAME% pam-scan-gui
