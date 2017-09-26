#!/usr/bin/env bash

echo ">>> runing the clipSearch to identify miRNA-target interactions "
../bin/clipSearch testGenome.fa testGenome.fa.fai testMir.fa testPeak.bed >test_clipSearch_mtis.txt 2>clipSearch.log.txt
