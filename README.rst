.. highlight:: tcsh

====================================================
Analysis code for natural scenes aperture experiment
====================================================

Preprocessing
=============

Subject-level
-------------

Filesystem
~~~~~~~~~~

1. Make the subject's directory structure::

    mkdir -p sXXXX/{analysis/,fmap,func/run{01,02,03,04,05,06,07,08,09,10,11,12},loc_analysis,logs,mvpa,reg}

2. Copy the subject's runtime logfiles to the ``logs`` directory.

3. Make symlinks named ``raw`` in each functional run directory (exp and loc) that link to the location of its associated raw DICOM directory::

    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/MR-such_and_such raw

4. Similarly, make symlinks named ``mag-raw`` and ``ph-raw`` in each fieldmap directory that link to the locations of the fieldmap acquisition::

    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/MR-SEyada mag-raw
    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/PH-SEyada ph-raw

5. Make a local copy of the AFNI/SUMA base anatomical that we can use for alignment::

    3dcopy \
       {$SUBJECTS_DIR}/{$SUBJ_ID}/SUMA/{$SUBJ_ID}_SurfVol+orig \
       reg/{$SUBJ_ID}_ns_aperture-anat+orig


Conversion
~~~~~~~~~~

Converts from the raw scanner format to a set of 4D NIFTI files::

    ns_aperture_preproc sXXXX convert

After execution, open up each NIFTI file and inspect for image quality and look at the summary image to see how much movement there was.


Fieldmap preparation
~~~~~~~~~~~~~~~~~~~~

Prepares the fieldmap::

    ns_aperture_preproc SXXXX fieldmap


Correction
~~~~~~~~~~
Applies a motion and distortion correction procedure::

    ns_aperture_preproc sXXXX mc_unwarp

After execution, open up the summary NIFTI file to check that most of the motion has been removed.
To verify that the unwarping has worked correctly:

* Run ``fslview``.
* Load the original or corrected image from a given run.
* Add the magnitude image from the fieldmap as an overlay.
* Notice the geometric distortions in the functional data.
* Add the undistorted image as an overlay, and hide the uncorrected image.
* Toggle the visibility of the undistorted image, and verify that the geometry now aligns well with that of the fieldmap's magnitude image.



Anatomical registration
~~~~~~~~~~~~~~~~~~~~~~~

First, make a copy of the mean functional::

    cd reg
    3dcopy ../func/sXXXX_ns_aperture-mean.nii sXXXX_ns_aperture-mean+orig

Now, we want to calculate some transformation parameters that will get the two images into rough register.
This will give the automated algorithm a good starting point.

* Start AFNI, from within the ``reg`` directory.
* Set the reference anatomical as the underlay.
* Position the crosshairs at a landmark on the brain. I like to use the most posterior portion of the occipital lobe, on the right side (in the image). Note down the three position values in the AFNI window (    in mm). Say they are ``[ 100, 50, 50 ]``.
* Then, change the underlay (or overlay, if you prefer) to the mean functional.
* Position the crosshairs at the same landmark as you used for the anatomical. The position might now be ``[ 20, 10, -20 ]``.
* Calculate ( reference anatomical positions - functional positions ), elementwise. In this example, that would give ``[ 80, 40, 70 ]``.
* Update the subject's configuration structure to include the estimate.

Then run::

    ns_aperture_preproc sXXXX sess_reg


Surface projection
~~~~~~~~~~~~~~~~~~
Projects the functional images to a standardised cortical surface, averaging between the white and pial surfaces::

    ns_aperture_preproc sXXXX vol_to_surf


Smoothing
~~~~~~~~~

Applies a small amount of spatial smoothing to the timecourses, calculated along a surface that is the average of the white and pial surfaces::

    ns_aperture_preproc sXXXX smooth


Temporal filtering
~~~~~~~~~~~~~~~~~~

High-pass filters each run's timecourse, for use in the MVPA procedure only::

    ns_aperture_preproc sXXXX filter_tc


Mask creation
~~~~~~~~~~~~~

Calculates a mask surface that specifies the nodes that are nonzero in each functional run::

    ns_aperture_preproc sXXXX surf_mask


Group-level
-----------

Filesystem
~~~~~~~~~~

Run::

    cd /labs/olmanlab/Data7T/NatSceneAperture
    mkdir -p group_data/{cluster_sim,loc,mvpa,ret}


Average anatomical
~~~~~~~~~~~~~~~~~~

Create an anatomical that is the average of all the subjects in the experiment::

    make_average_subject \
      -subjects s1000 s1008 s1011 s1021 s1032 \
      -out ns_aperture_avg \
      -sd-out /labs/olmanlab/NatSceneAperture/group_data/


Mask creation
~~~~~~~~~~~~~

Create a mask that indicates the nodes that are nonzero in all the individual subject masks::

    ns_aperture_group_analysis group_mask



Cluster simulation
~~~~~~~~~~~~~~~~~~

Runs a Monte-Carlo simulation to determine the FWE cluster threshold for each hemisphere::

    ns_aperture_group_analysis clust_sim

After running the above, edit the two scripts to replace ``Surf_A`` in the ``SurfClust`` call with ``midway``. Then (takes a while)::

    script_lh.sh
    script_rh.sh

Retinotopy
~~~~~~~~~~

Average each subject's wedge and ring maps, in a phase-sensitive way::

    ns_aperture_group_analysis ret_std


Searchlight preparation
~~~~~~~~~~~~~~~~~~~~~~~

Identifies the nodes associated with a searchlight around each node. The single-subject ``mvpa_prep`` needs to have occured, for the representative subject, prior to running this. Takes ages::

    ns_aperture_group_analysis mvpa_node_prep


Univariate experiment analysis
==============================


Subject-level
-------------

Design preparation
~~~~~~~~~~~~~~~~~~
Computes the experimental design from the logfiles::

    ns_aperture_analysis sXXXX design_prep


GLM
~~~

Estimate the GLM::

    ns_aperture_analysis sXXXX glm


Cluster summary
~~~~~~~~~~~~~~~

After the group-level clustering has been done, run::

    ns_aperture_analysis sXXXX coh_clust_summ


Group-level
-----------

Height threshold
~~~~~~~~~~~~~~~~

Runs a one-sample t-test on the subject beta weights::

    ns_aperture_group_analysis coh_test


Cluster threshold
~~~~~~~~~~~~~~~~~

To apply the cluster threshold::

    ns_aperture_group_analysis coh_clust


Cluster summary
~~~~~~~~~~~~~~~

To print out a summary of the cluster beta statistics::

    ns_aperture_group_analysis coh_effect_size



Univariate localiser analysis
=============================


Subject-level
-------------

Design preparation
~~~~~~~~~~~~~~~~~~

Generate the design info::

    ns_aperture_analysis sXXXX loc_design_prep

GLM
~~~

Execute the GLM::

    ns_aperture_analysis sXXXX loc_glm


Group-level
-----------

Height threshold
~~~~~~~~~~~~~~~~

Runs a one-sample t-test on the ( either > 0 ) regressor::

    ns_aperture_analysis loc_test


Cluster threshold
~~~~~~~~~~~~~~~~~


Multivariate analysis
=====================

Subject-level
-------------

Design and data preparation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This saves the node info, the condition info, and the z-scored block data for a given subject. Run::

    ns_aperture_analysis sXXXX mvpa_prep


Classification
~~~~~~~~~~~~~~

Group-level
-----------

Height threshold
~~~~~~~~~~~~~~~~

Cluster threshold
~~~~~~~~~~~~~~~~~

