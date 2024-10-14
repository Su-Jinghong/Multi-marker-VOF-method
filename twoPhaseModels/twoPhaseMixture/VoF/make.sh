#!/bin/sh

chmod 777 -R ./;
dir=`pwd`

for((i=4;i<=16;i++));
do
   cd $dir
   pwd
   cp -rf alphaEqn2.H alphaEqn$i.H;
   sed -i s:bubble2:bubble$i:g alphaEqn$i.H;
done

