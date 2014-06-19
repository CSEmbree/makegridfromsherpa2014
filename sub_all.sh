#!/bin/bash

export STEER=$1

export DIR=`pwd`

if (test $STEER) then
 echo steering is $STEER number of events $2
else
 echo no steering found
 exit
fi

echo -------> launch phase space generation
echo
rm $DIR/outputtop/*
#./makegridfromsherpa $1 -*_*  $2
./makegridfromsherpa $1 -B
./makegridfromsherpa $1 -R
echo -------> phase space generation done
#
#
echo -------> launch grid generation
echo
#cp -r $DIR/outputtopsave/* $DIR/outputtop
# 
./makegridfromsherpa $1 -B | tee $1-B.log
./makegridfromsherpa $1 -R | tee $1-R.log
#
#./makegridfromsherpa $1 -B_*  $2 > $1-B_*.log
#./makegridfromsherpa $1 -R_*  $2 > $1-R_*.log
#
#./makegridfromsherpa $1 -B*  $2 > $1-B_*.log
#./makegridfromsherpa $1 -R*  $2 > $1-R_*.log
echo -------> launch grid generation

#echo
#./makegridfromsherpa $1 $2