@echo off
setlocal
set ROOT=%~dp0
set PY=C:\software\Python311\python.exe
set PPM=
set OUT=

for /f "delims=" %%f in ('dir /b /s /o-d "%ROOT%..\out\binary.ppm" 2^>nul') do (
    set PPM=%%f
    goto :found
)
:found

if "%PPM%"=="" (
    echo [x] binary.ppm not found.
    echo     Build and run Assignment5 first.
    pause
    exit /b 1
)

if exist "%ROOT%result.png" del /f /q "%ROOT%result.png" >nul 2>&1

"%PY%" "%ROOT%tools\ppm2png.py" "%PPM%" "%ROOT%result.png"
if errorlevel 1 (
    echo [x] convert failed
    pause
    exit /b 1
)

for /f "delims=" %%f in ('dir /b /s /o-d "%ROOT%result*.png" 2^>nul') do (
    set OUT=%%f
    goto :opened
)
:opened

start "" "%OUT%"
