#!/bin/bash -eu

# SPDX-FileCopyrightText: 2025 Uwe Fechner
#
# SPDX-License-Identifier: MIT

rm -rf /tmp/testdepot
mkdir /tmp/testdepot
rm -rf /tmp/test
mkdir /tmp/test
cd /tmp/test
export JULIA_DEPOT_PATH=/tmp/testdepot 
julia --project="." -e "using Pkg; pkg\"add KiteModels@0.6.3\"; using KiteModels"
cd ..
rm -rf /tmp/test2   
mkdir /tmp/test2
cd /tmp/test2
julia --project="." -e "using Pkg; pkg\"add KiteModels\"; using KiteModels; pkg\"status\""
cd ..

