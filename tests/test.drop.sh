#!/usr/bin/env bash

pip install "${HOME}/Dropbox_p/repos/nwkit"
wd="${HOME}/Dropbox_p/repos/nwkit/tests/"
cd ${wd}

echo "Testing internal node info drop from nhx"
nwkit drop --infile ../data/generax.nhx --outfile - --target intnode --support yes --name yes

echo "Testing file to stdout"
nwkit drop --infile ../data/OG0001999.dated.pruned.nwk --outfile - --target root --length yes
nwkit drop --infile ../data/OG0001999.dated.pruned.nwk --outfile - --target root --length yes --fill 0

echo "Testing stdin to stdout"
cat ../data/OG0001999.dated.pruned.nwk \
| nwkit drop --target root --length yes
