ECHO ON

pushd %~dp0..
"%~dp0muscle5.1.win64.exe" -disperse align\ensemble.efa -log align\disperse.log
popd

PAUSE
