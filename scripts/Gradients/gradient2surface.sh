# Get schaefer200-16subcort atlas in surface space (for spin-test)

# Set up FreeSurfer environment
export FREESURFER_HOME=/Users/au704655/Documents/Packages/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

# covert atlas to surface space also
mri_vol2surf --src /Users/au704655/Documents/Body_wandering/Data/gradients/atlas_schaefer200_7net_3d.nii.gz --out schaefer200subcort_atlas.L.mgh --hemi lh --regheader fsaverage
mri_vol2surf --src /Users/au704655/Documents/Body_wandering/Data/gradients/atlas_schaefer200_7net_3d.nii.gz --out schaefer200subcort_atlas.R.mgh --hemi rh --regheader fsaverage

