#!/bin/bash
extra_docker_opts="--shm-size=$(($(stat -f --format="%b*%s" /dev/shm/)))"
docker run $extra_docker_opts -i --tty=false -a STDIN -a STDOUT -a STDERR --user $EUID:$EUID \
       --volume=$(pwd):$(pwd) --workdir=$(pwd) $@
