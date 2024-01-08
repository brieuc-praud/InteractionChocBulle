#!/usr/bin/env bash

# =====
# Use the isentropic Euler vortex test case
# to perform a convergence analysis
# =====

# Parameter file for the program
PARAMS=parameters.dat
# Output to this file
ERRFILE=error.dat
# $1 is the name of the executable
EXE="$1"

# Starting point for convergence analysis
REFNx=50
REFNy=50
REFCFL=0.25
PWRbasis=1.1
REFtmax=10.
# Number of points for the convergence analysis
NMAX=4

# Get the line of each parameter
NNx=$(cat -n $PARAMS | grep "Nx" | cut -f 1)
NNy=$(cat -n $PARAMS | grep "Ny" | cut -f 1)
NCFL=$(cat -n $PARAMS | grep "CFL" | cut -f 1)
Ntmax=$(cat -n $PARAMS | grep "tmax" | cut -f 1)
NOutputMod=$(cat -n $PARAMS | grep "output_modulo" | cut -f 1)
# Macro to change the parameters of interest
changeParameters()
{
    awk "NR == $NNx { print \"Nx\", $1 }\
         NR == $NNy { print \"Ny\", $2 }\
         NR == $NCFL { print \"CFL\", $3 }\
         NR == $Ntmax { print \"tmax\", $4 }\
         NR == $NOutputMod { print \"output_modulo\", $5 }\
         NR != $NNx && NR != $NNy && NR != $NCFL && NR != $Ntmax && NR != $NOutputMod { print \$0 }" $PARAMS > tmp.txt
    mv tmp.txt $PARAMS
}

# Store the old parameters to restore them after the analysis
cp $PARAMS $PARAMS.tmp
function cleanup {
    mv $PARAMS.tmp $PARAMS
}
trap cleanup EXIT # restore the config file when exiting
# Create a new file (erase the content if it already exists)
> $ERRFILE
for i in $(seq 1 $NMAX)
do
    echo -ne "  $i/$NMAX\r"

    Nx=$(echo "$REFNx * $PWRbasis^$i" | bc -l)
    Nx=${Nx%.*} #Conversion to an integer
    Ny=$(echo "$REFNy * $PWRbasis^$i" | bc -l)
    Ny=${Ny%.*}
    CFL=$(echo "$REFCFL / $PWRbasis^$i" | bc -l)
    OutputMod=-1

    dx=$(echo "1./$Nx" | bc -l)

    changeParameters $Nx $Ny $CFL $REFtmax $OutputMod

    output=$(exec "./$EXE")
    err=$(echo "$output" | tail -n 1 | cut -d':' -f 2 | xargs)
    echo "$dx $err" >> $ERRFILE
done
