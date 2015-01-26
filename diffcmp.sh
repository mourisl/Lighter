#!/bin/sh

for f in `ls ../lighter_git/`
do
	echo $f
	diff $f ../lighter_git
done
