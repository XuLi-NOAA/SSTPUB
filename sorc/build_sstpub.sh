#! /usr/bin/env bash
set -eux

module purge
source ../modulefiles/sstpub.wcoss_dell_p3  
module list

cd ./sstpub.fd
./makefile.sh
