ECHO ON

cd %~dp0
copy "..\paml\clean_aligned_paml.phylip" "clean_aligned_paml.phylip"
copy "..\paml\radseq_pruned_unrooted.tre" "radseq_pruned_unrooted.tre"
"%~dp0codeml.exe"

PAUSE
