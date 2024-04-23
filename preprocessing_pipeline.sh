#!/bin/bash
# Preprocessing steps to 

## Path 
# data dir
DATADIR=/local_path/

#Synb0-DISCO-local
SYNB0_DIR=/local_path/Synb0-DISCO/local
export PATH=$PATH:${SYNB0_DIR}

# CATO
CATODIR=/local_path/CATO/3.2.2-macos
export PATH=$PATH:${CATODIR}
MCRDIR=/local_path/MATLAB_Runtime/v93

cd $DATADIR

### run SynB0-DISCO for subjects without reverse encoded b0 scan
for sub in `cat` synb0_subjects.txt; do
    echo ""
    echo "-----------------------"
    echo "Processing: $sub   "     
    # get dataset
    dataset=$(grep -o -E "driv|vici|wrap" <<< "$sub")

    # set paths
    T1=$DATADIR/data-${dataset}/rawdata/${sub}/ses-01/anat/${sub}_ses-01_T1w.nii.gz
    B0=$DATADIR/data-${dataset}/derivatives/synb0-DISCO/ses-01/${sub}/${sub}_b0.nii.gz
    OUT=$DATADIR/data-${dataset}/derivatives/synb0-DISCO/ses-01/${sub}
    DWI=$DATADIR/data-${SET}/rawdata/${SUB}/ses-01/dwi/${SUB}_ses-01_dwi.nii.gz

    # topup requires even dimensions (add zero slice if odd)
    check_nii_dims.sh $B0 $B0
    check_nii_dims.sh $DWI $DWI

    # run synb0
    synb0-disco_local.sh -t $T1 -b $B0 -o $OUT -i

done


## run CATO for all datasets
# - every dataset has its own config file
# - datasets without reverse-phase encoding epi images have their own acq file
# - datasets with reverse-phase encoding epi images use the default CATO acq file

# loop over all datasets
for SET in driv vici wrap ercp; do
    # set freesurfer subjects dir to correct path
    SUBJECTS_DIR=/localpath/data-${SET}/derivatives/freesurfer/ses-01
    
    # loop over subjects in the dataset
    for SUB in `grep $SET subjects.txt`; do

        # subject path
        SUBPATH=$DATADIR/data-${SET}/rawdata/${SUB}       

        # run cato
        structural_pipeline --configurationFile=$DATADIR/conf_${SET}.json \
        -s $SUBPATH -m $MCRDIR
        
    done
done