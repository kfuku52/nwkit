#!/usr/bin/env bash

pip install '/Users/kf/Dropbox_w/repos/nwkit'

wd="/Users/kf/Dropbox_w/repos/nwkit/tests/"
cd ${wd}


nwkit drop --infile ../data/OG0001999.dated.pruned.nwk --outfile - --target root --length yes

