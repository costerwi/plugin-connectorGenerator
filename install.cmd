@echo off

set plugin=connectorGenerator

if defined HOME (
    set abaqus_plugins=%HOME%\abaqus_plugins
) else (
    set abaqus_plugins=%HOMEDRIVE%%HOMEPATH%\abaqus_plugins
)

echo This script copies the Abaqus CAE plugin "%plugin%" into "%abaqus_plugins%"

if not exist "%abaqus_plugins%" mkdir "%abaqus_plugins%"

set destination=%abaqus_plugins%\%plugin%\
if "%~dp0"=="%destination%" (
    echo ERROR: Already installed.
    pause
    exit 1
)

if not exist "%destination%" mkdir "%destination%"
copy /Y "%~dp0\*.*" "%destination%"
if ERRORLEVEL 0 (
    echo Success! Restart Abaqus CAE and check Plugin-ins menu.
) else (
    echo Something went wrong.
    echo Create an issue on https://github.com/costerwi
)
pause
