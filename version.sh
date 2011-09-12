#!/bin/zsh
# ccs, to be called by cmake to get the 
# current build version
# useful useful useful
v=`git describe --tags --long`
echo -n $v