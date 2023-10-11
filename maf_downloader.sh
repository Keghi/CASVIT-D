#!/bin/bash

# set the URL to the website
url="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz470way/maf/"

# use curl to get a list of all files on the website
files=$(curl -s $url | grep -o '<a href=['"'"'"][^"'"'"']*gz['"'"'"]' | sed -e 's/^<a href=["'"'"']//' -e 's/["'"'"']$//')

# loop through each file and download it using wget
for file in $files; do
  wget -b $url$file
done

for file in *.gz
do
  md5sum "$file" >> md5sum_output.txt
  gunzip "$file"
done
