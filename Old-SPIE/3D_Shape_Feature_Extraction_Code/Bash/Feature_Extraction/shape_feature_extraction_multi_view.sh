#!/bin/bash

#STEPS (copy and paste these into command window)
#1) login to rider.case.edu (using Putty)
#
#2) set path:
#		 cd /mnt/projects/CSE_BME_AXM788/home/cxa220/Projects/Rectal_Radiology/Code/Bash/Feature_Extraction_Modified/
#2) convert scripts from windows to unix:
#		 dos2unix shape_feature_extraction_CCF.sh
#3) call main script!
#		 ./shape_feature_extraction_CCF.sh
#4) check status of jobs
#		 squeue -u <username>


#DATA=/mnt/projects/CSE_BME_AXM788/home/cxa220/Projects/Rectal_Radiology/Data/CCFRectalCA
DATA="/scratch/pbsjobs/job.17255576.hpc/Coronal_Resampled"
#PROG=/mnt/Data7/jta35/Code/ShapeFeatures3D_HPC/build/ShapeFeatures3D
PROG="/home/tgd15/Shape_Features/Shape_Feature_Extraction_Code/C/ShapeFeatures3D_HPC/build/ShapeFeatures3D"
#RESULTS=/mnt/projects/CSE_BME_AXM788/home/cxa220/Projects/Rectal_Radiology/Data/CCFRectalCA
RESULTS="/scratch/pbsjobs/job.17255576.hpc/Coronal_Resampled"

#check patient number
cutoff_1=10
cutoff_2=100

#labels
declare -a arr=("lumen" "rw" "fat")

## a loop through the above array
for roi in "${arr[@]}"
do
	#echo "$i"
	for ((j=1; j<=180; j++)) # move patient by patient
	do

	    if [[ "$j" -lt "$cutoff_1" ]];then #if j<10, we need the "0" in front
		    patient=UH-RectalCA\-\0\0$j
	    fi

	    if [[ "$j" -ge "$cutoff_1" ]] && [[ "$j" -lt "$cutoff_2" ]];then
		    patient=UH-RectalCA\-\0$j
	    fi

	    if [[ "$j" -ge "$cutoff_2" ]];then
		    patient=UH-RectalCA\-$j
	    fi

	    #extract shape from isotropic volume!
	    vol=$DATA/$patient/vol/vol\.mha
	    mask=$DATA/$patient/masks/mask_$roi\.mha

	    if [[ -e $mask ]]
	    then
		# Call to Extract 3D Features
		# <$PATH/ShapeFeatures3D> <LabelImage> <IntensityImage> <OutputFeatureFile.txt>
				$PROG $mask $vol "$RESULTS/$patient/shape_feats_label_$roi.txt"
	    else
				echo "Path does not exist $mask"
		fi
	done
done
