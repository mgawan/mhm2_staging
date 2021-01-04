#!/bin/bash

set -e


BBTOOLS_IMG=${BBTOOLS_IMG:=bryce911/bbtools:latest}
BBTOOLS=
isdocker=0
if [ -x "$(which shifter)" ]
then
  shifterimg lookup $BBTOOLS_IMG || shifterimg pull $BBTOOLS_IMG
  BBTOOLS="shifter --image=$BBTOOLS_IMG"
elif [ -x "$(which docker)" ]
then
  isdocker=1
  docker pull $BBTOOLS_IMG
  BBTOOLS="docker run $BBTOOLS_IMG"
else
  echo "Please install docker or shifter"
  exit 1
fi

USAGE="$0 <reads.fastq>[.gz] [ read2.fastq[.gz] ]

This script will perform some basic trimming on the input fastq to remove
illumina sequencing artifacts and adaptors.

If the data is paired in two files, the resulting file will be interleaved
as a single uncompressed file, output to the current directory as : *-bbqc.fq

You may want to set TMPDIR to a place with enough disk to manipulate your FASTQ files

"

tmpdir=$(mktemp -d)
tmpdir2=$tmpdir

input=$1
input2=$2
if [ ! -f "$input" ]
then
  echo "$USAGE"
  echo "You must specify an input file"
  echo
  exit 1
fi

in2="int"
if [ -f "$input2" ]
then
  in2="in2=$input2"
fi

if [ "$isdocker" == "1" ]
then
  BBTOOLS="docker run -i --tty=false -a STDIN -a STDOUT -a STDERR --user $(id -u):$(id -g) -w $(pwd) --volume=$(pwd):$(pwd)  --mount=type=bind,source=$tmpdir,target=$tmpdir ${BBTOOLS_IMG}"
fi

shortin=${input##*/}
outbase=${shortin%.gz}
outbase=${outbase%.fastq}
outbase=${outbase%.fq}
bbdukopts=${BBDUK_OPTS:=}

catorzcat=zcat
if [ "${input%.gz}" == "${input}" ]
then
  catorzcat=cat
fi

LFS=${LFS:=lfs}
LFS=$(which ${LFS} || /bin/true)
set_dir_striping()
{
  local stripe=$1
  local dir=$2
  echo "set_dir_striping ${stripe} ${dir}"
  if [ -x "${LFS}" ]
  then
     lfsdf=$(${LFS} df ${dir})
     if [ -z "${lfsdf}" ]
     then
         echo "LUSTRE is not on ${dir}"
     else
       if [ "$stripe" == "72" ]
       then
          stripe=$(( $(echo "${lfsdf}" | grep -c _UUID) * 8 / 10))
          echo "using stripe ${stripe} instead of 72"
       fi

       if [ ${stripe} -le 1 ]
       then
           stripe=1
       fi
       echo "Setting lustre stripe to ${stripe} on ${dir}"
       ${LFS} setstripe -c ${stripe} ${dir} || true
     fi
  else
     echo "LFS not detected, skipping stripe count $1 on $2"
  fi
}

if [ ! -f ${outbase}-bbqc.fq.gz ]
then
  tmpoutdir=${outbase}.tmp
  rm -rf ${tmpoutdir}
  mkdir -p ${tmpoutdir}
  set_dir_striping 72 ${tmpoutdir}
  set_dir_striping 1 ${tmpdir}

  #Trim adapters
  chastity=$(${catorzcat} $input | head -1 | grep ' [12]:[YN]:.:' || true)
  if [ -n "${chastity}" ] ; then chastity=chastityfilter ; fi

  out=${tmpoutdir}/${outbase##*/}-bbqc.fq

  $BBTOOLS bash -c "bbduk.sh -Xmx1g in=$input $in2 out=stdout.fq ow qin=auto qout=33 ktrim=r k=23 mink=11 hdist=1 ref=/bbmap/resources/adapters_no_transposase.fa.gz tbo tpe minlen=35 $chastity maxns=1 ftm=5 $bbdukopts \| bbduk.sh -Xmx1g in=stdin.fq int out=$out int ow k=27 hdist=1 qtrim=f ref=/bbmap/resources/sequencing_artifacts.fa.gz,/bbmap/resources/phix174_ill.ref.fa.gz minlen=35"
#  $BBTOOLS bash -c "bbduk.sh -Xmx1g in=$input $in2 out=stdout.fq ow qin=auto qout=33 ktrim=r k=23 mink=11 hdist=1 ref=adapters tbo tpe minlen=35 $chastity maxns=1 ftm=5 $bbdukopts \| bbduk.sh -Xmx1g in=stdin.fq int out=$out int ow k=27 hdist=1 qtrim=f ref=/bbmap/resources/sequencing_artifacts.fa.gz,/bbmap/resources/phix174_ill.ref.fa.gz minlen=35"

  mv ${tmpoutdir}/* .
  rmdir ${tmpoutdir}

fi

rm -rf ${tmpdir}

echo "

Finished.  Your trimmed fastq files should be here:

$(ls -l ${outbase}-bbqc.fq)

"

