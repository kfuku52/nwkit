#!/usr/bin/env bash

pip install '/Users/kf/Dropbox_w/repos/nwkit'

wd="/Users/kf/Dropbox_w/repos/nwkit/tests/"
cd ${wd}

echo "Testing file to stdout"
nwkit drop --infile ../data/OG0001999.dated.pruned.nwk --outfile - --target root --length yes
nwkit drop --infile ../data/OG0001999.dated.pruned.nwk --outfile - --target root --length yes --fill 0


echo "Testing stdin to stdout"
cat ../data/OG0001999.dated.pruned.nwk \
| nwkit drop --target root --length yes
