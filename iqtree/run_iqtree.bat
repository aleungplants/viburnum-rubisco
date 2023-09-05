ECHO ON

pushd %~dp0
"%~dp0iqtree2.exe" -s "clean_aligned.phylip" -m MFP -B 5000 -T 3 -redo
"%~dp0codeml.exe"
popd

PAUSE
