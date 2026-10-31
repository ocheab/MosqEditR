@echo off
setlocal
set "URL=https://ngs.sanger.ac.uk/production/ag1000g/phase2/AR1/variation/main/vcf/pass/ag1000g.phase2.ar1.pass.2R.vcf.gz"
echo Checking that Windows Subsystem for Linux is installed...
wsl.exe --status
if errorlevel 1 echo WSL is not ready. See README_START_HERE.txt.
echo.
echo Checking small HTTP byte-range request with HEAD only. NO chromosome download.
curl.exe -sS -I -H "Range: bytes=0-31" "%URL%"
echo.
echo Look for HTTP 206 Partial Content and Content-Range, ideally. If server responds 200, remote random access may not work.
pause
