# mhmxx #

mhmxx is a UPC++ version of [MetaHipMer](https://sites.google.com/lbl.gov/exabiome/downloads?authuser=0).

This code relies on UPC++, which can be obtained from

https://bitbucket.org/berkeleylab/upcxx/wiki/Home


## Building

To build, simply run

`./build.sh Release`

or

`./build.sh Debug`

Once it has been built once, you can rebuild with

`./build.sh`

and it will build using the previously chosen setting (Release or Debug)

You can also run

`./build.sh clean`

to really start from scratch.

## Running

A typical command line to run (e.g. on 10 nodes each with 24 processors) from the root source directory is:

`upcxx-run -shared-heap 5% -n 240 -- ./build/mhmxx --cores-per-node 24 -r <READS.fastq> -i <insert_size_avg:insert_size_stddev> -k 21,33,55,77,99 -s 99,33`

This will run the application in the current directory. We suggest that you create a separate output directory and run from there.

Run with `-h` to see options.

## Cori notes:

To build and run on [Cori](https://docs.nersc.gov/systems/cori/), you'll need the following modules:

`module load cmake`  

`module load upcxx`

If building for KNL, make sure to first do

`module switch craype-haswell craype-mic-knl`

Also, currently this does not build for the Intel compiler on Cori, so use Cray (e.g. or GNU

For Cori, it is recommended to stripe all files on the Lustre file system to ensure adequate IO performance, e.g.

`lfs setstripe -c -1 reads.fastq`
