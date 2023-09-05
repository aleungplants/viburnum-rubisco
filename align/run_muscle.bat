ECHO ON

pushd %~dp0..
"%~dp0muscle5.1.win64.exe" -align download\clean.fasta -output align\clean_aligned.fasta -threads 7
python "%~dp0fasta_to_phylip.py"
"%~dp0muscle5.1.win64.exe" -align download\clean.fasta -stratified -output align\ensemble.efa -threads 7
"%~dp0muscle5.1.win64.exe" -disperse align\ensemble.efa -log align\disperse.log
popd

PAUSE
