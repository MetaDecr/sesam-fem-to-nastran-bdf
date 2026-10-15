@echo off
rem ============================================================
rem  fem2bdf GUI launcher - starts the local server and opens
rem  the browser. Keep this file ASCII-only with CRLF endings.
rem  errorlevel 9009 = interpreter not found; only then fall
rem  back to the next Python command.
rem ============================================================
cd /d %~dp0

python server.py --open
if errorlevel 9009 py -3 server.py --open
if errorlevel 9009 python3 server.py --open
if errorlevel 9009 (
  echo.
  echo [ERROR] No working Python found ^(tried: python, py -3, python3^).
)
pause
