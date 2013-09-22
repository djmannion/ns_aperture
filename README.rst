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

Surface projection
~~~~~~~~~~~~~~~~~~
Projects the functional images to the cortical surface::

    ns_aperture_preproc sXXXX vol_to_surf


Smoothing
~~~~~~~~~

Temporal filtering
~~~~~~~~~~~~~~~~~~

Mask creation
~~~~~~~~~~~~~

Group-level
-----------

Mask creation
~~~~~~~~~~~~~

Cluster simulation
~~~~~~~~~~~~~~~~~~

Retinotopy
~~~~~~~~~~

Searchlight preparation
~~~~~~~~~~~~~~~~~~~~~~~


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

Group-level
-----------

Height threshold
~~~~~~~~~~~~~~~~

Cluster threshold
~~~~~~~~~~~~~~~~~


Univariate localiser analysis
=============================


Subject-level
-------------

Design preparation
~~~~~~~~~~~~~~~~~~

GLM
~~~

Group-level
-----------

Height threshold
~~~~~~~~~~~~~~~~

Cluster threshold
~~~~~~~~~~~~~~~~~


Multivariate analysis
=====================

Subject-level
-------------

Design and data preparation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Classification
~~~~~~~~~~~~~~

Group-level
-----------

Height threshold
~~~~~~~~~~~~~~~~

Cluster threshold
~~~~~~~~~~~~~~~~~

