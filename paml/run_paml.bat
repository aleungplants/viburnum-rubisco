ECHO ON

cd %~dp0
"C:\Program Files\R\R-4.2.2\bin\Rscript.exe" unroot_tree.R
"%~dp0codeml.exe"

PAUSE
