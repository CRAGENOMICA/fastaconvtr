@test "fa2fa.ex01 should give same results as input" {
  run ../bin/fastaconvtr -F fasta -f fasta -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_01.fa -n ./chr10.txt
  [ "$status" -eq 0 ]
  [ "$output" = "" ]
  [ -f "./100Kchr10_fa2fa_01.fa" ]
}