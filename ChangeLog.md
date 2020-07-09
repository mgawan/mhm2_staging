# ChangeLog.md for MetaHipMer version 2 (aka mhmxx)


This is the ChangeLog for MetaHipMer with development at [bitbucket](https://bitbucket.org/berkeleylab/mhmxx)

### 0.1.2 2020-07-06
   * Semi-stable alpha release
   * Added skeleton for CI building on cori and hulk development server
   * Fixed re-build issues with upcxx-utils
   * Added warning on build and install when submodule is out of date
   * Various fixes for Summit build including vector support for SSW alignment
   * Added support for more cluster schedulers
   * Added --paired-reads support and modified Fastq class to support two-file libraries
   * Fixed LUSTRE striping on per_thread rank logs


### 0.1.1 Before 2020-07-01
   * Complete rewrite of MetaHipMer 1.0 using only UPC++
