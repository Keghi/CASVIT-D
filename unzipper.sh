#!/bin/bash

for file in *.gz
do
  md5sum "$file" >> md5sum_output.txt
  gunzip "$file"
done