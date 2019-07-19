#!/bin/bash

fn="run_cable_spatial.py"
of="temp.py"

sed 's/standard/hydraulics/g' $fn > $of
mv $of $fn

sed -e 's#../../src/trunk/trunk/#../../src/trunk_DESICA_PFTs/trunk_DESICA_PFTs/#g' $fn > $of
mv $of $fn
