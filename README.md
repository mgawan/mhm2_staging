# mhmxx #

mhmxx is a UPC++ version of [MetaHipMer](https://sites.google.com/lbl.gov/exabiome/downloads?authuser=0).

This code relies on UPC++, which can be obtained from

https://bitbucket.org/berkeleylab/upcxx/wiki/Home


## Building

To build, first set your environment:

`export MHMXX_BUILD_ENV=.env-cori-knl/gnu.sh`

and then simply run:

`./build.sh Release`

or:

`./build.sh Debug`

The script will install the binaries by default into the `install/bin` subdirectory in the repository root directory. To set a different install 
directory, set the environment variable `MHMXX_INSTALL_PATH`, e.g.:

`MHMXX_INSTALL_PATH=$SCRATCH/mhmxx-install ./build.sh Release`

Once mhmxx has been built once, you can rebuild with

`./build.sh`

and it will build using the previously chosen setting (Release or Debug)

You can also run

`./build.sh clean`

to really start from scratch.

You can also use cmake directly:

`mkdir -p build && cd build && cmake -DCMAKE_INSTALL_PREFIX=path-to-install .. && make -j all install`

You'll need to first set the environment, e.g.:

`source .env-cori-knl/gnu.sh`


## Running


A typical command line to run is:

`mhmxx.py -r <READS_lib1.fastq>,<READS_lib2.fastq>`

This will create a new output directory that contains the results. Note that mhmxx requires interleaved paired reads.

Run with `-h` to see the various options.

## Cori notes:

To build and run on [Cori](https://docs.nersc.gov/systems/cori/), you'll need the upcxx module and the cmake module.

All six permuations of gnu, cray and intel environments on haswell and knl hardware are supported
by sourcing one of these environments:

`.env-cori-haswell/cray.sh  .env-cori-haswell/gnu.sh  .env-cori-haswell/intel.sh  .env-cori-knl/cray.sh  .env-cori-knl/gnu.sh  .env-cori-knl/intel.sh`

To use these with the build.sh script, simply run as, e.g.:

`MHMXX_BUILD_ENV=.env-cori/env-knl-gnu.sh ./build.sh Release`

It is recommended to use either PrgEnv-gnu or PrgEnv-cray. Builds with the Intel compiler run very slowly and currently there seems to be some bug that causes the execution to hang. This is not present with the other compilers.

It is recommended to stripe all input files on the Lustre file system to ensure adequate I/O performance, e.g.

`lfs setstripe -c -1 reads.fastq`

The output directory (created using the `-o` parameter) is automatically striped.
