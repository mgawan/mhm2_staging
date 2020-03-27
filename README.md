# mhmxx #

mhmxx is a UPC++ version of [MetaHipMer](https://sites.google.com/lbl.gov/exabiome/downloads?authuser=0).

This code relies on UPC++, which can be obtained from

https://bitbucket.org/berkeleylab/upcxx/wiki/Home


## Building

To build, simply run

`./build.sh Release`

or

`./build.sh Debug`

These will install the binaries by default into the `bin` subdirectory in the repository root directory. To set a different install 
directory, set the environment variable `MHMXX_INSTALL_PATH`, e.g.:

`MHMXX_INSTALL_PATH=$SCRATCH/mhmxx-install ./build.sh Release`

Once mhmxx has been built once, you can rebuild with

`./build.sh`

and it will build using the previously chosen setting (Release or Debug)

You can also run

`./build.sh clean`

to really start from scratch.

## Running


A typical command line to run is:

`mhmxx.py -r <READS_lib1.fastq>,<READS_lib2.fastq> -i <insert_size_avg:insert_size_stddev>`

This will create a new output directory that contains the results.

Run with `-h` to see the various options.

## Cori notes:

To build and run on [Cori](https://docs.nersc.gov/systems/cori/), you'll need the upcxx module:

`module load upcxx`

If building for KNL, make sure to first do

`module switch craype-haswell craype-mic-knl`

For Cori, it is recommended to stripe all input files on the Lustre file system to ensure adequate I/O performance, e.g.

`lfs setstripe -c -1 reads.fastq`

The output directory (created using the `-o` parameter) is automatically striped.
