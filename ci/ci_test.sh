#!/bin/bash

mhmxx_dir=$1
data_dir=$2

rm -rf /dev/shm/test-as
$mhmxx_dir/install/bin/mhmxx.py -r $data_dir/arctic_sample_0.fq -o /dev/shm/test-as --checkpoint=no
$mhmxx_dir/ci/check_asm_quality.py --asm-dir /dev/shm/test-as --expected-quals $mhmxx_dir/ci/good-arctic-sample0.txt --refs $data_dir/arcticsynth-refs.fa 
