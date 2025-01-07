#To compile:
#gcc ./sources/*.c -lm -o ./bin/fastaconvtr -Wall -O3 -g -O0 -lz
sh ./build.sh
cp ./build/fastaconvtr ./bin

cd ./Examples

../bin/fastaconvtr -h > ../fastaconvtr_help.txt

#FASTA TO TFASTA 
echo --------------------------------------------------------------------------------------------------
echo fasta to tfasta: Useful to run mstatspop in sliding windows mode. Also generates a weighting file.
echo --------------------------------------------------------------------------------------------------

echo
echo fa2tfa.ex01
echo ../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10.tfa.gz -n ./chr10.txt
../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10.tfa.gz -n ./chr10.txt
echo
echo fa2tfa.ex01b masking several regions with Ns
echo ../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_01B.tfa.gz -m ./coord_100Kb.txt -n ./chr10.txt
../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_01B.tfa.gz -m ./coord_100Kb.txt -n ./chr10.txt
echo
echo fa2tfa.ex02
echo ../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_02.tfa.gz  -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_02.tfa.gz  -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -n ./chr10.txt
echo
echo fa2tfa.ex03 should give same results than previous
echo ../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_03.tfa.gz  -E ./100Kchr10_fa2tfa_02.tfa.gz_npops1_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_03.tfa.gz  -E ./100Kchr10_fa2tfa_02.tfa.gz_npops1_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
echo
echo fa2tfa.ex02B
echo ../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_02B.tfa.gz -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n ./chr10.txt  -n ./chr10.txt
../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_02B.tfa.gz -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n ./chr10.txt  -n ./chr10.txt
echo
echo fa2tfa.ex03B should give same results than previous
echo ../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_03B.tfa.gz -E ./100Kchr10_fa2tfa_02B.tfa.gz_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_03B.tfa.gz -E ./100Kchr10_fa2tfa_02B.tfa.gz_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
echo
echo fa2tfa.ex04
echo ../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_04.tfa.gz -p 2 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_04.tfa.gz -p 2 -n ./chr10.txt
echo
echo fa2tfa.ex05
echo ../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_05.tfa.gz -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -p 2 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_05.tfa.gz -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -p 2 -n ./chr10.txt
echo
echo fa2tfa.ex06 should give same results than previous
echo ../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_06.tfa.gz -E ./100Kchr10_fa2tfa_05.tfa.gz_npops1_nsam84_nonsynonymous_max_IncludeMissingVariantsmhits_NOoutg_ploidy2_TFAv2_WEIGHTS.gz -p 2 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_06.tfa.gz -E ./100Kchr10_fa2tfa_05.tfa.gz_npops1_nsam84_nonsynonymous_max_IncludeMissingVariantsmhits_NOoutg_ploidy2_TFAv2_WEIGHTS.gz -p 2 -n ./chr10.txt
echo
echo fa2tfa.ex07
echo ../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_07.tfa.gz -N 2 40 2 -G 1 -u 1 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_07.tfa.gz -N 2 40 2 -G 1 -u 1 -n ./chr10.txt
echo
echo fa2tfa.ex08
echo ../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_08.tfa.gz -N 2 40 2 -G 1 -u 1 -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -n ./chr10.txt
../bin/fastaconvtr -F fasta -f tfasta  -i ./100Kchr10.fa -o ./100Kchr10_fa2tfa_08.tfa.gz -N 2 40 2 -G 1 -u 1 -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -n ./chr10.txt
echo

#FASTA TO FASTA
echo --------------------------------------------------------------------------------------------------
echo fasta to fasta: Useful for concatenate different regions from coordenates file. Also generating weighting file from GFF file
echo --------------------------------------------------------------------------------------------------
echo
echo fa2fa.ex01 should give same results than input
echo ../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_01.fa  -n ./chr10.txt
../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_01.fa -n ./chr10.txt
echo
echo fa2fa.ex01B should give same results than input
echo ../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10_fa2fa_01.fa -o ./100Kchr10_fa2fa_01B.fa  -n ./chr10.txt
../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10_fa2fa_01.fa -o ./100Kchr10_fa2fa_01B.fa  -n ./chr10.txt
echo
echo fa2fa.ex02
echo ../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_02.fa -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_02.fa -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -n ./chr10.txt
echo
echo fa2fa.ex02B
echo ../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_02B.fa -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_02B.fa -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n ./chr10.txt
echo
echo fa2fa.ex03 should give same results than previous
echo ../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_03.fa -E ./100Kchr10_fa2fa_02.fa_npops1_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_03.fa -E ./100Kchr10_fa2fa_02.fa_npops1_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
echo
echo fa2fa.ex03B should give same results than previous
echo ../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_03B.fa -E ./100Kchr10_fa2fa_02B.fa_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_03B.fa -E ./100Kchr10_fa2fa_02B.fa_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
echo
echo fa2fa.ex04
echo ../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_04.fa -W ./coord_100Kb.txt -n ./chr10.txt
../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_04.fa -W ./coord_100Kb.txt -n ./chr10.txt
echo
echo fa2fa.ex05
echo ../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_05.fa -p 2 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_05.fa -p 2 -n ./chr10.txt
echo
echo fa2fa.ex06
echo ../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_06.fa -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -p 2 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_06.fa -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -p 2 -n ./chr10.txt
echo
echo fa2fa.ex07 should give same results than previous
echo ../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_07.fa -E ./100Kchr10_fa2fa_06.fa_npops1_nsam84_nonsynonymous_max_IncludeMissingVariantsmhits_NOoutg_ploidy2_TFAv2_WEIGHTS.gz -p 2 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f fasta   -i ./100Kchr10.fa -o ./100Kchr10_fa2fa_07.fa -E ./100Kchr10_fa2fa_06.fa_npops1_nsam84_nonsynonymous_max_IncludeMissingVariantsmhits_NOoutg_ploidy2_TFAv2_WEIGHTS.gz -p 2 -n ./chr10.txt
echo

#FASTA TO MS
echo --------------------------------------------------------------------------------------------------
echo fasta to ms: Useful to generate mask files in case doing simulations. msfile contains only variants in full-non-missing positions
echo --------------------------------------------------------------------------------------------------
echo
echo fa2ms.ex01
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_01.ms.txt -n ./chr10.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_01.ms.txt -n ./chr10.txt
echo
echo fa2ms.ex02
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_02.ms.txt -N 1 42 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_02.ms.txt -N 1 42 -n ./chr10.txt
echo
echo fa2ms.ex03
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_03.ms.txt -N 2 40 2 -G 1 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_03.ms.txt -N 2 40 2 -G 1 -n ./chr10.txt
echo
echo fa2ms.ex03B
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_03B.ms.txt -N 2 80 4 -G 1 -p 2 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_03B.ms.txt -N 2 80 4 -G 1 -p 2 -n ./chr10.txt
echo
echo fa2ms.ex04
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_04.ms.txt -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n ./chr10.txt > ./100Kchr10_fa2ms_04.ms.txt.log.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_04.ms.txt -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n ./chr10.txt > ./100Kchr10_fa2ms_04.ms.txt.log.txt
echo
echo fa2ms.ex04b
echo ../bin/fastaconvtr -F fasta -f ms  -N 2 40 2 -G 1 -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_04b.ms.txt -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n ./chr10.txt > ./100Kchr10_fa2ms_04b.ms.txt.log.txt
../bin/fastaconvtr -F fasta -f ms  -N 2 40 2 -G 1 -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_04b.ms.txt -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n ./chr10.txt > ./100Kchr10_fa2ms_04b.ms.txt.log.txt
echo
echo
echo fa2ms.ex05 should give same results than previous
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_05.ms.txt -E ./100Kchr10_fa2ms_04.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_05.ms.txt -E ./100Kchr10_fa2ms_04.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
echo
echo fa2ms.ex06
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_06.ms.txt -W ./coord_100Kb.txt -n ./chr10.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_06.ms.txt -W ./coord_100Kb.txt -n ./chr10.txt
echo
echo fa2ms.ex06b should give same results than previous
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_06b.ms.txt -w 10000 -s 20000 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_06b.ms.txt -w 10000 -s 20000 -n ./chr10.txt
echo
echo fa2ms.ex07
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_07.ms.txt -N 2 40 2 -w 10000 -s 20000 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_07.ms.txt -N 2 40 2 -w 10000 -s 20000 -n ./chr10.txt
echo
echo fa2ms.ex08
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_08.ms.txt -P 0 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_08.ms.txt -P 0 -n ./chr10.txt
echo
echo fa2ms.ex09
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_09.ms.txt -p 2 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_09.ms.txt -p 2 -n ./chr10.txt
echo
echo fa2ms.ex10
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_10.ms.txt -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -p 2 -n ./chr10.txt > ./100Kchr10_fa2ms_10.ms.txt.log.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_10.ms.txt -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -p 2 -n ./chr10.txt > ./100Kchr10_fa2ms_10.ms.txt.log.txt
echo
echo fa2ms.ex11 should give same results than previous
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_11.ms.txt -E ./100Kchr10_fa2ms_10.ms.txt_npops1_nsam84_nonsynonymous_max_ExcludeMissingmhits_NOoutg_ploidy2_TFAv2_WEIGHTS.gz -p 2 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_11.ms.txt -E ./100Kchr10_fa2ms_10.ms.txt_npops1_nsam84_nonsynonymous_max_ExcludeMissingmhits_NOoutg_ploidy2_TFAv2_WEIGHTS.gz -p 2 -n ./chr10.txt
echo
echo fa2ms.ex12
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_12.ms.txt -W ./coord_100Kb.txt -p 2 -n ./chr10.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_12.ms.txt -W ./coord_100Kb.txt -p 2 -n ./chr10.txt
echo
echo fa2ms.ex13
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_13.ms.txt -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -W ./coord_100Kb.txt -P 0 -n ./chr10.txt > ./100Kchr10_fa2ms_13.ms.txt.log.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_13.ms.txt -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -W ./coord_100Kb.txt -P 0 -n ./chr10.txt > ./100Kchr10_fa2ms_13.ms.txt.log.txt
echo
echo fa2ms.ex14 should give same results than previous
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_14.ms.txt -E ./100Kchr10_fa2ms_13.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -P 0 -p 2 -W ./coord_100Kb.txt -n ./chr10.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_14.ms.txt -E ./100Kchr10_fa2ms_13.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -P 0 -p 2 -W ./coord_100Kb.txt -n ./chr10.txt
echo
echo fa2ms.ex15
echo ../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_15.ms.txt -N 3 40 40 4 -G 1 -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -p 2 -n ./chr10.txt > ./100Kchr10_fa2ms_15.ms.txt.log.txt
../bin/fastaconvtr -F fasta -f ms      -i ./100Kchr10.fa -o ./100Kchr10_fa2ms_15.ms.txt -N 3 40 40 4 -G 1 -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -p 2 -n ./chr10.txt > ./100Kchr10_fa2ms_15.ms.txt.log.txt
echo

#TFASTA TO MS
echo --------------------------------------------------------------------------------------------------
echo tfasta to ms: Useful to generate mask files in case doing simulations. msfile contains only variants in full-non-missing positions
echo --------------------------------------------------------------------------------------------------
echo
echo tfa2ms.ex01
echo ../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_01.ms.txt -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_01.ms.txt -n ./chr10.txt
echo
echo tfa2ms.ex02
echo ../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_02.ms.txt -N 1 42 -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_02.ms.txt -N 1 42 -n ./chr10.txt
echo
echo tfa2ms.ex03
echo ../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_03.ms.txt -N 2 40 2 -G 1 -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_03.ms.txt -N 2 40 2 -G 1 -n ./chr10.txt
echo
echo tfa2ms.ex04
echo ../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_04.ms.txt -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n ./chr10.txt > ./100Kchr10_tfa2ms_04.ms.txt.log.txt
../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_04.ms.txt -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n ./chr10.txt > ./100Kchr10_tfa2ms_04.ms.txt.log.txt
echo
echo tfa2ms.ex05 should give same results than previous
echo ../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_05.ms.txt -E ./100Kchr10_tfa2ms_04.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingmhits_NOoutg_ploidy1_WEIGHTS.txt -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_05.ms.txt -E ./100Kchr10_tfa2ms_04.ms.txt_npops1_nsam42_nonsynonymous_max_ExcludeMissingmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
echo
echo tfa2ms.ex06
echo ../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_06.ms.txt -W ./coord_100Kb.txt -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_06.ms.txt -W ./coord_100Kb.txt -n ./chr10.txt
echo
echo tfa2ms.ex07 should give same results than previous
echo ../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_07.ms.txt -w 10000 -s 20000 -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_07.ms.txt -w 10000 -s 20000 -n ./chr10.txt
echo
echo tfa2ms.ex08
echo ../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_08.ms.txt -P 0 -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_08.ms.txt -P 0 -n ./chr10.txt
echo
echo tfa2ms.ex09
echo ../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_09.ms.txt -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -W ./coord_100Kb.txt -P 0 -n ./chr10.txt > ./100Kchr10_tfa2ms_09.ms.txt.log.txt
../bin/fastaconvtr -F tfasta -f ms     -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2ms_09.ms.txt -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -W ./coord_100Kb.txt -P 0 -n ./chr10.txt > ./100Kchr10_tfa2ms_09.ms.txt.log.txt
echo
echo tfa2ms.ex10: checking multiple scaffolds
echo ../bin/fastaconvtr -F tfasta -f ms     -i ./100Kallchr.tfa.bgz -o ./100Kchr10_tfa2ms_10.ms.txt -P 0 -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f ms     -i ./100Kallchr.tfa.bgz -o ./100Kchr10_tfa2ms_10.ms.txt -P 0 -n ./chr10.txt
echo
echo tfa2ms.ex11: checking multiple scaffolds coordinates file
echo ../bin/fastaconvtr -F tfasta -f ms     -i ./100Kallchr.tfa.bgz -o ./100Kchr10_tfa2ms_11.ms.txt -W ./coord_100Kb_allchr.txt -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f ms     -i ./100Kallchr.tfa.bgz -o ./100Kchr10_tfa2ms_11.ms.txt -W ./coord_100Kb_allchr.txt -n ./chr10.txt
echo

#TFASTA TO FASTA
echo --------------------------------------------------------------------------------------------------
echo tfasta to fasta: Useful to generate a weighting file from GFF file
echo --------------------------------------------------------------------------------------------------
echo
echo tfa2fa.ex01
echo ../bin/fastaconvtr -F tfasta -f fasta  -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2fa_01.fa -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f fasta  -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2fa_01.fa -n ./chr10.txt
echo
echo tfa2fa.ex01 masking several regions with Ns
echo ../bin/fastaconvtr -F tfasta -f fasta  -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2fa_01B.fa -m ./coord_100Kb.txt -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f fasta  -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2fa_01B.fa -m ./coord_100Kb.txt -n ./chr10.txt
echo
echo tfa2fa.ex02
echo ../bin/fastaconvtr -F tfasta -f fasta  -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2fa_02.fa -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -n ./chr10.txt > ./100Kchr10_tfa2fa_02.fa.log.txt
../bin/fastaconvtr -F tfasta -f fasta  -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2fa_02.fa -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -n ./chr10.txt > ./100Kchr10_tfa2fa_02.fa.log.txt
echo
echo tfa2fa.ex02B
echo ../bin/fastaconvtr -F tfasta -f fasta  -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2fa_02B.fa -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n ./chr10.txt > ./100Kchr10_tfa2fa_02B.fa.log.txt
../bin/fastaconvtr -F tfasta -f fasta  -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2fa_02B.fa -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n ./chr10.txt > ./100Kchr10_tfa2fa_02B.fa.log.txt
echo
echo tfa2fa.ex03 should give same results than previous
echo ../bin/fastaconvtr -F tfasta -f fasta  -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2fa_03.fa -E ./100Kchr10_tfa2fa_02.fa_npops1_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f fasta  -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2fa_03.fa -E ./100Kchr10_tfa2fa_02.fa_npops1_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
echo
echo tfa2fa.ex03B should give same results than previous
echo ../bin/fastaconvtr -F tfasta -f fasta  -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2fa_03B.fa -E ./100Kchr10_tfa2fa_02B.fa_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f fasta  -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2fa_03B.fa -E ./100Kchr10_tfa2fa_02B.fa_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
echo
echo include multiple scaffolds tfasta: check simple, include gtf, coord, mask...
echo
echo tfa2tfa.ex04 using multiple chromosomes and output multiple chromosomes
echo ../bin/fastaconvtr -F tfasta -f fasta -i ./100Kallchr.tfa.bgz -o ./100Kchr1012_tfa2fa_04.fa -n ./chr1012.txt
../bin/fastaconvtr -F tfasta -f fasta -i ./100Kallchr.tfa.bgz -o ./100Kchr1012_tfa2fa_04.fa -n ./chr1012.txt
echo
echo tfa2tfa.ex05 using as input multiple chromosomes + gtf
echo ../bin/fastaconvtr -F tfasta -f fasta -i ./100Kallchr.tfa.bgz -o ./100Kchr12_tfa2fa_05.fa -n ./chr12.txt -g ./100Kallchr.gtf nonsynonymous Nuclear_Universal -c max -u 1
../bin/fastaconvtr -F tfasta -f fasta -i ./100Kallchr.tfa.bgz -o ./100Kchr12_tfa2fa_05.fa -n ./chr12.txt -g ./100Kallchr.gtf nonsynonymous Nuclear_Universal -c max -u 1
echo
echo tfa2tfa.ex06 using as input multiple chromosomes + coord
echo ../bin/fastaconvtr -F tfasta -f fasta -i ./100Kallchr.tfa.bgz -o ./100Kchr12_tfa2fa_06.fa -n ./chr12.txt -W ./coord_100Kb_allchr.txt -u 1
../bin/fastaconvtr -F tfasta -f fasta -i ./100Kallchr.tfa.bgz -o ./100Kchr12_tfa2fa_06.fa -n ./chr12.txt -W ./coord_100Kb_allchr.txt -u 1
echo
echo tfa2tfa.ex07 using as input multiple chromosomes + mask
echo ../bin/fastaconvtr -F tfasta -f fasta -i ./100Kallchr.tfa.bgz -o ./100Kchr12_tfa2fa_07.fa -n ./chr12.txt -m ./coord_100Kb_allchr.txt -u 1
../bin/fastaconvtr -F tfasta -f fasta -i ./100Kallchr.tfa.bgz -o ./100Kchr12_tfa2fa_07.fa -n ./chr12.txt -m ./coord_100Kb_allchr.txt -u 1
echo

#TFASTA TO TFASTA
echo --------------------------------------------------------------------------------------------------
echo tfasta to tfasta: Useful for generating weighting file from GFF file
echo --------------------------------------------------------------------------------------------------
echo
echo tfa2tfa.ex02
echo ../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2tfa_02.tfa.gz -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -n ./chr10.txt > ./100Kchr10_tfa2tfa_02.tfa.gz.log.txt
../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2tfa_02.tfa.gz -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -n ./chr10.txt > ./100Kchr10_tfa2tfa_02.tfa.gz.log.txt
echo
echo tfa2tfa.ex02B
echo ../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2tfa_02B.tfa.gz -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n ./chr10.txt > ./100Kchr10_tfa2tfa_02B.tfa.gz.log.txt
../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2tfa_02B.tfa.gz -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 0 -n ./chr10.txt > ./100Kchr10_tfa2tfa_02B.tfa.gz.log.txt
echo
echo tfa2tfa.ex03 should give same results than previous
echo ../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2tfa_03.tfa.gz -E ./100Kchr10_tfa2tfa_02.tfa.gz_npops1_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2tfa_03.tfa.gz -E ./100Kchr10_tfa2tfa_02.tfa.gz_npops1_nsam42_nonsynonymous_max_IncludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
echo
echo tfa2tfa.ex03B should give same results than previous
echo ../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2tfa_03B.tfa.gz -E ./100Kchr10_tfa2tfa_02B.tfa.gz_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2tfa_03B.tfa.gz -E ./100Kchr10_tfa2tfa_02B.tfa.gz_npops1_nsam42_nonsynonymous_max_ExcludeMissingVariantsmhits_NOoutg_ploidy1_TFAv2_WEIGHTS.gz -n ./chr10.txt
echo
echo tfa2tfa.ex03c
echo ../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2tfa_03c.tfa.gz -N 3 40 40 4 -G 1 -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -p 1 -n ./chr10.txt 
../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kchr10.tfa.gz -o ./100Kchr10_tfa2tfa_03c.tfa.gz -N 3 40 40 4 -G 1 -g ./100Kchr10.gtf nonsynonymous Nuclear_Universal -c max -u 1 -p 1 -n ./chr10.txt
echo
echo tfa2tfa.ex04 using multiple chromosomes
echo ../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kallchr.tfa.gz -o ./100Kallchr.tfa.bgz  -n ./chr101214.txt
../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kallchr.tfa.gz -o ./100Kallchr.tfa.bgz  -n ./chr101214.txt
echo
echo tfa2tfa.ex05 using multiple chromosomes
echo ../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kallchr.tfa.gz -o ./100Kallchr10i14.tfa.gz  -n ./chr1014.txt
../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kallchr.tfa.gz -o ./100Kallchr10i14.tfa.gz  -n ./chr1014.txt
echo
echo tfa2tfa.ex06 using multiple chromosomes
echo ../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kallchr.tfa.gz -o ./100Kallchr12i14.tfa.gz  -n ./chr1214.txt
../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kallchr.tfa.gz -o ./100Kallchr12i14.tfa.gz  -n ./chr1214.txt
echo
echo tfa2tfa.ex07 using multiple chromosomes
echo ../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kallchr.tfa.gz -o ./100Kallchr07.tfa.gz  -g ./100Kallchr.gtf nonsynonymous Nuclear_Universal -c max -u 1 -n ./chr101214.txt -N 2 40 2 -G 1
../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kallchr.tfa.gz -o ./100Kallchr07.tfa.gz -g ./100Kallchr.gtf nonsynonymous Nuclear_Universal -c max -u 1 -n ./chr101214.txt -N 2 40 2 -G 1
echo
echo tfa2tfa.ex08 using multiple chromosomes and masking multiple chromosomes
echo ../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kallchr.tfa.gz -o ./100Kallchr08.tfa.gz  -n ./chr101214.txt -m ./coord_100Kb_allchr.txt
../bin/fastaconvtr -F tfasta -f tfasta -i ./100Kallchr.tfa.gz -o ./100Kallchr08.tfa.gz  -n ./chr101214.txt -m ./coord_100Kb_allchr.txt 
echo 
