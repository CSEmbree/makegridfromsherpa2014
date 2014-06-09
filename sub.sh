#!/bin/bash

export STEER=$1

export DIR=`pwd`

if (test $STEER) then
 echo steering is $STEER number of events $2
else
 echo no steering found
 exit
fi

echo -------> assuming phase space generation is done
#
#
echo -------> launch grid generation
echo
#cp -r $DIR/outputtopsave/* $DIR/outputtop
# 
./makegridfromsherpa $1 -B  $2 > $1-B.log
./makegridfromsherpa $1 -R  $2 > $1-R.log
#
#./makegridfromsherpa $1 -B_*  $2 > $1-B_*.log
#./makegridfromsherpa $1 -R_*  $2 > $1-R_*.log
#
#./makegridfromsherpa $1 -B*  $2 > $1-B_*.log
#./makegridfromsherpa $1 -R*  $2 > $1-R_*.log
echo -------> launch grid generation

#echo
#./makegridfromsherpa $1 $2