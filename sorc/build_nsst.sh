#! /usr/bin/env bash
set -eux

module purge
source ../modulefiles/nsst.wcoss_dell_p3  
module list

mkdir -p ../exec

cd ./nsst.fd
./makefile.sh
