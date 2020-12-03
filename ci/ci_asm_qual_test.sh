#!/bin/bash

mhm2_install_dir=$(dirname $(dirname $(realpath $0) ) )
if [ -n "$1" ] || [ ! -x ${mhm2_install_dir}/bin/ci_asm_qual_test.sh ]
then
  mhm2_install_dir=$(realpath $1)
fi

refs=arcticsynth-refs.fa
if [ ! -f "$refs" ]; then
    wget https://portal.nersc.gov/project/hipmer/MetaHipMer_datasets_12_2019/ArcticSynth/samples/${refs}.gz
    gunzip ${refs}.gz
fi
reads=arctic_sample_0.fq
if [ ! -f "$reads" ]; then
    wget https://portal.nersc.gov/project/hipmer/MetaHipMer_datasets_12_2019/ArcticSynth/samples/${reads}.gz
    gunzip ${reads}.gz
fi 

wd=`pwd`
test_dir=$wd/test-arctic-sample0
rm -rf $test_dir
$mhm2_install_dir/bin/mhm2.py -r $reads -o $test_dir --checkpoint=no --post-asm-align --post-asm-abd
status=$?
if [ $status -ne 0 ]
then
  echo "MHM2 failed! - $status"
  exit $status
fi
$mhm2_install_dir/bin/check_asm_quality.py --asm-dir $test_dir --expected-quals $mhm2_install_dir/share/good-arctic-sample0.txt --refs $wd/$refs
status=$?
if [ $status -ne 0 ]
then
  echo "check_asm_quality.py failed! - $status"
  exit $status
fi
