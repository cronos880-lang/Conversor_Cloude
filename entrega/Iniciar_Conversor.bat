@echo off
setlocal
cd /d "%~dp0"
title Conversor 3D STL a STEP
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0servidor_local.ps1"
if errorlevel 1 pause
