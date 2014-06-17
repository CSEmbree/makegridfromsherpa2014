#!/bin/bash

export STEER=$1
export DIR=`pwd`
#
#
#
if (test $STEER) then
 echo steering is $STEER number of events $2
else
 echo no steering found
 exit
fi
#
#
#
echo -------> launch phase space generation
echo
#rm $DIR/outputtop/*
#
./makegridfromsherpa $1 -B | tee $1-B-phase.log
./makegridfromsherpa $1 -R | tee $1-R-phase.log
#
echo -------> finished phase space generation
#
#
#
echo -------> launch grid generation
echo
#cp $DIR/outputtop/* $DIR/outputtop-save-temp/
# 
./makegridfromsherpa $1 -B | tee $1-B.log
./makegridfromsherpa $1 -R | tee $1-R.log
#
echo -------> finshed grid generation
echo