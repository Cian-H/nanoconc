#!/bin/bash
# This bash script will pull and compile the other implementations of the
# bhmie algorithm (found at http://scatterlib.wikidot.com/mie)
# These implementations will then be tested and benchmarked against the local
# implementation to ensure that they are correct and (ideally) faster.

# First, let's create a base URL for the scatterlib wiki codebases
base_url="http://scatterlib.wikidot.com/local--files/codes/"

# Next, let's create a list of the codebases we want to pull
codebases=("bhmie-f.zip" "bhmie-c.zip")

# if bhmie_dir containers bhmie-c/bhmie.so, bhmie-f/bhmie.so, and bhmie-f/bhmie_f77.so then we can skip the build
bhmie_dir=$1
if [ -f $bhmie_dir/bhmie-c/bhmie.so ] && [ -f $bhmie_dir/bhmie-f/bhmie.so ] && [ -f $bhmie_dir/bhmie-f/bhmie_f77.so ]; then
    echo "bhmie-c/bhmie.so, bhmie-f/bhmie.so, and bhmie-f/bhmie_f77.so already exist. Skipping build."
    exit 0
fi

# Now, let's pull the codebases
for codebase in ${codebases[@]}; do
    wget -O $bhmie_dir/$codebase $base_url/$codebase
    if [[ $codebase == *.zip ]]; then
        unzip $bhmie_dir/$codebase -d $bhmie_dir/$codebase
    fi
done

# Then, let's extract any .zip files to a directory of the same name (without the .zip, obviously)
for codebase in ${codebases[@]}; do
    if [[ $codebase == *.zip ]]; then
        unzip $bhmie_dir/$codebase -d $bhmie_dir/${codebase%.zip}
    fi
done

# If the folder already existed in the archive, move that folder one level up
# and remove the now empty folder
for codebase in ${codebases[@]}; do
    if [[ $codebase == *.zip ]]; then
        folder=${codebase%.zip}
        if [ -d $bhmie_dir/$folder/$folder ]; then
            mv $bhmie_dir/$folder/$folder/* $bhmie_dir/$folder
            rmdir $bhmie_dir/$folder/$folder
        fi
    fi
done

# Next, if this has all succeeded we can delete the .zip files
for codebase in ${codebases[@]}; do
    if [[ $codebase == *.zip ]]; then
        rm $bhmie_dir/$codebase
    fi
done

# And, finally, we can compile the C, and Fortran implementations
cd $bhmie_dir
gcc -g -shared -fPIC -o bhmie-c/bhmie.so bhmie-c/bhmie.c bhmie-c/complex.c bhmie-c/nrutil.c -lm -Wno-builtin-declaration-mismatch -Wno-implicit-function-declaration
gfortran -g -shared -fPIC -o bhmie-f/bhmie.so bhmie-f/bhmie.f
gfortran -g -shared -fPIC -o bhmie-f/bhmie_f77.so bhmie-f/bhmie_f77.f