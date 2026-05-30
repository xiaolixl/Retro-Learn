@echo off
chcp 65001 >nul 2>&1
echo ============================================
echo  SimpRetro - Dependency Installation
echo ============================================
echo.

where python >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python not found. Please install Python 3.11+ and add it to PATH.
    pause
    exit /b 1
)

echo [1/6] Upgrading pip...
python -m pip install --upgrade pip
echo.

echo [2/6] Installing numpy...
python -m pip install numpy==1.26.4
echo.

echo [3/6] Installing matplotlib...
python -m pip install matplotlib==3.8.3
echo.

echo [4/6] Installing tqdm...
python -m pip install tqdm==4.67.1
echo.

echo [5/6] Installing rdkit-pypi and rdchiral...
python -m pip install rdkit-pypi==2024.9.5
if errorlevel 1 (
    echo rdkit-pypi version pin failed, trying latest...
    python -m pip install rdkit-pypi
)
python -m pip install rdchiral==1.1.0
echo.

echo [6/6] Verifying all packages...
python -c "import rdkit; import rdchiral; import numpy; import matplotlib; import tqdm; print('All packages imported successfully!')"
echo.

if errorlevel 1 (
    echo.
    echo [WARNING] Some packages failed to import. Check the errors above.
) else (
    echo.
    echo ============================================
    echo  Installation Complete!
    echo ============================================
)
pause
