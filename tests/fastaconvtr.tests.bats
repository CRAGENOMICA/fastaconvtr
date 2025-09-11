load ./base.tests.bats
load '../node_modules/bats-support/load'
load '../node_modules/bats-assert/load'
# bats test_tags=tag:bin
@test "fastaconvtr ok" {
  run fastaconvtr -h
  assert_success
  [ "${lines[2]}" = "#fastaconvtr " ] 
  #"#fastaconvtr v0.1beta (20220907) Sebastian E. Ramos-Onsins." ]
}


# FASTA TO TFASTA tests with tag fa2tfa
# echo --------------------------------------------------------------------------------------------------
# echo fasta to tfasta: Useful to run mstatspop in sliding windows mode. Also generates a weighting file.
# echo --------------------------------------------------------------------------------------------------

# bats test_tags=tag:fa2tfa
@test "fa2tfa ex01" {
  ## rm $TEST_OUTPUT/*
  run fastaconvtr -F fasta -f tfasta  -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10.tfa.gz -n $TEST_FILES_DIR/chr10.txt
  assert_success
  #run gzip -df $TEST_OUTPUT/100Kchr10.tfa.gz 
  # diff skiping first line
  #run diff -q <(tail -n +2 $TEST_FILES_DIR/100Kchr10.tfa) <(tail -n +2 $TEST_OUTPUT/100Kchr10.tfa)
  # run diff -q $TEST_FILES_DIR/100Kchr10.tfa $TEST_OUTPUT/100Kchr10.tfa 
  #assert_success
}

# bats test_tags=tag:fa2tfa
@test "fa2tfa ex01b masking several regions with Ns" {
  run fastaconvtr -F fasta -f tfasta  -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2tfa_01B.tfa.gz -m $TEST_FILES_DIR/coord_100Kb.txt -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2tfa
@test "fa2tfa ex02" {
  # ...
  run fastaconvtr -F fasta -f tfasta  -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2tfa_02.tfa.gz  -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success


}

# bats test_tags=tag:fa2tfa
@test "fa2tfa ex03" {
  # ...
  run fastaconvtr -F fasta -f tfasta -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2tfa_03B.tfa.gz -E $TEST_FILES_DIR/100Kchr10_fa2tfa_02B.NonSyn.WEIGHTS.gz -n $TEST_FILES_DIR/chr10.txt
  assert_success


}

# bats test_tags=tag:fa2tfa
@test "fa2tfa ex02B" {
  # ...
  run fastaconvtr -F fasta -f tfasta -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2tfa_02B.tfa.gz -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2tfa
@test "fa2tfa ex03B should give same results of fa2tfa ex02B" {
  # ...
  run fastaconvtr  -F fasta -f tfasta -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2tfa_03B.tfa.gz -E $TEST_FILES_DIR/100Kchr10_fa2tfa_02B.NonSyn.WEIGHTS.gz -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2tfa
@test "fa2tfa ex04" {
  # ...
  run fastaconvtr  -F fasta -f tfasta -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2tfa_04.tfa.gz -p 2 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2tfa
@test "fa2tfa ex05" {
  # ...
  run fastaconvtr -F fasta -f tfasta -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2tfa_05.tfa.gz -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -p 2 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2tfa
@test "fa2tfa ex06 should give same results of fa2tfa ex05" {
  # ...
  run fastaconvtr  -F fasta -f tfasta -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2tfa_06.tfa.gz -E $TEST_FILES_DIR/100Kchr10_fa2tfa_05.NonSyn.WEIGHTS.gz -p 2 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2tfa
@test "fa2tfa ex07" {
  # ...
  run fastaconvtr -F fasta -f tfasta -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2tfa_07.tfa.gz -N 2 40 2 -G 1 -u 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}


# bats test_tags=tag:fa2tfa
@test "fa2tfa ex08" {
  # ...
  run fastaconvtr -F fasta -f tfasta -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2tfa_08.tfa.gz -N 2 40 2 -G 1 -u 1 -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -n $TEST_FILES_DIR/chr10.txt
  assert_success

}


#FASTA TO FASTA
# --------------------------------------------------------------------------------------------------
# fasta to fasta: Useful for concatenate different regions from coordenates file. Also generating weighting file from GFF file
# --------------------------------------------------------------------------------------------------
# bats test_tags=tag:fa2fa
@test "fa2fa ex01" {
  # echo fa2fa.ex01 should give same results than input
  run fastaconvtr -F fasta -f fasta  -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2fa_01.fa -n $TEST_FILES_DIR/chr10.txt
  assert_success

}


# bats test_tags=tag:fa2fa
@test "fa2fa ex01B" {
  # echo fa2fa.ex01B should give same results than input
  run fastaconvtr  -F fasta -f fasta   -i $TEST_FILES_DIR/100Kchr10_fa2fa_01.fa -o $TEST_OUTPUT/100Kchr10_fa2fa_01B.fa  -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2fa ex02" {
  run fastaconvtr  -F fasta -f fasta  -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2fa_02.fa -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2fa ex02B" {
  # ...
  run fastaconvtr -F fasta -f fasta  -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2fa_02B.fa -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2fa ex03 same results as fa2fa ex02B" {
  # ...
  run fastaconvtr  -F fasta -f fasta  -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2fa_03.fa -E $TEST_FILES_DIR/100Kchr10_fa2fa_02.NonSyn.WEIGHTS.gz -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2fa ex03B same results as fa2fa ex02B" {
  # ...
  run fastaconvtr -F fasta -f fasta  -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2fa_03B.fa -E $TEST_FILES_DIR/100Kchr10_fa2fa_02B.NonSyn.WEIGHTS.gz -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2fa ex04" {
  # ...
  run fastaconvtr -F fasta -f fasta  -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2fa_04.fa -W $TEST_FILES_DIR/coord_100Kb.txt -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2fa ex05" {
  # ...
  run fastaconvtr  -F fasta -f fasta  -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2fa_05.fa -p 2 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2fa ex06" {
  # ...
  run fastaconvtr -F fasta -f fasta  -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2fa_06.fa -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -p 2 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}


# bats test_tags=tag:fa2fa
@test "fa2fa ex07 as results of fa2fa ex06" {
  # ...
  run fastaconvtr -F fasta -f fasta  -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2fa_07.fa -E $TEST_FILES_DIR/100Kchr10_fa2fa_06.NonSyn.WEIGHTS.gz -p 2 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}


#FASTA TO MS
# --------------------------------------------------------------------------------------------------
# fasta to ms: Useful to generate mask files in case doing simulations. msfile contains only variants in full-non-missing positions
# --------------------------------------------------------------------------------------------------


# bats test_tags=tag:fa2fa
@test "fa2ms ex01" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_01.ms.txt -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex02" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_02.ms.txt -N 1 42 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex03" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_03.ms.txt -N 2 40 2 -G 1 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex03B" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_03B.ms.txt -N 2 80 4 -G 1 -p 2 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex04" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_04.ms.txt -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n $TEST_FILES_DIR/chr10.txt > $TEST_FILES_DIR/100Kchr10_fa2ms_04.ms.txt.log.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex04b" {
  # ...
  run fastaconvtr -F fasta -f ms  -N 2 40 2 -G 1 -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_04b.ms.txt -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n $TEST_FILES_DIR/chr10.txt > $TEST_FILES_DIR/100Kchr10_fa2ms_04b.ms.txt.log.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex05 same results as fa2ms ex04b" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_05.ms.txt -E $TEST_FILES_DIR/100Kchr10_fa2ms_04.NonSyn.WEIGHTS.gz -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex06" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_06.ms.txt -W $TEST_FILES_DIR/coord_100Kb.txt -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex06b same results as fa2ms ex06" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_06b.ms.txt -w 10000 -s 20000 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex07" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_07.ms.txt -N 2 40 2 -w 10000 -s 20000 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex08" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_08.ms.txt -P 0 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex09" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_09.ms.txt -p 2 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex10" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_10.ms.txt -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -p 2 -n $TEST_FILES_DIR/chr10.txt > $TEST_FILES_DIR/100Kchr10_fa2ms_10.ms.txt.log.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex11 same results as fa2ms ex10" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_11.ms.txt -E $TEST_FILES_DIR/100Kchr10_fa2ms_10.NonSyn.WEIGHTS.gz -p 2 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex12" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_12.ms.txt -W $TEST_FILES_DIR/coord_100Kb.txt -p 2 -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex13" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_13.ms.txt -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -W $TEST_FILES_DIR/coord_100Kb.txt -P 0 -n $TEST_FILES_DIR/chr10.txt > $TEST_FILES_DIR/100Kchr10_fa2ms_13.ms.txt.log.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex14 as results of fa2ms ex13" {
  # ...
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_14.ms.txt -E $TEST_FILES_DIR/100Kchr10_fa2ms_13.NonSyn.WEIGHTS.gz -P 0 -p 2 -W $TEST_FILES_DIR/coord_100Kb.txt -n $TEST_FILES_DIR/chr10.txt
  assert_success

}

# bats test_tags=tag:fa2fa
@test "fa2ms ex15" {
  # ... 
  run fastaconvtr -F fasta -f ms -i $TEST_FILES_DIR/100Kchr10.fa -o $TEST_OUTPUT/100Kchr10_fa2ms_15.ms.txt -N 3 40 40 4 -G 1 -g $TEST_FILES_DIR/100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -p 2 -n $TEST_FILES_DIR/chr10.txt > $TEST_FILES_DIR/100Kchr10_fa2ms_15.ms.txt.log.txt
  assert_success

}


# @test 'refute_output' {
#   run echo 'want'
#   refute_output 'want'
#   echo 'want' | assert_output
# }