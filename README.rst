====================================================
Analysis code for Glass pattern coherence experiment
====================================================

Requirements
============

- Python >= 2.7
- numpy
- matplotlib
- scipy
- fmri_tools (`http://www.bitbucket.org/djmannion/fmri_tools <http://www.bitbucket.org/djmannion/fmri_tools/>`_ )


Processing stages
=================

Prepare the filesystem
----------------------

1. Make the subject's directory structure::

    mkdir -p sXXXX/{anat,analysis,fmap/f{1,2},func/run{1,2,3,4},log,roi}

2. Copy the subject's runtime logfile to the ``log`` directory.

3. Make symlinks named ``raw`` in each functional run directory that link to the location of its associated raw DICOM directory::

    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/MR-such_and_such raw

4. Similarly, make symlinks named ``mag-raw`` and ``ph-raw`` in each fieldmap directory that link to the locations of the fieldmap acquisition::

    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/MR-SEyada mag-raw
    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/PH-SEyada ph-raw

5. Copy (not symlink) the subject's main high-res anatomical (skull stripped) from the main repository to the ``anat`` directorym using FSL::

    fslmaths /labs/olmanlab/Anatomy/sXXXX/sXXXX_stripped sXXXX_anat

  N.B. We want to use single NIFTI files, so before running the above you may need to run::

      setenv FSLOUTPUTTYPE NIFTI

6. Copy the relevant ROI MAT files from the visual localisers repository (the Gray view) to the ``roi`` directory.


Update the experiment information file
--------------------------------------

Edit ``get_subj_conf`` within ``glass_coherence/config.py`` and add the new subject's information.

For example::

    sXXXX = { "subj_id" : "sXXXX",
              "acq_date" : "YYYYMMDD",
              "n_runs" : 4,
              "n_fmaps" : 1,
              "comments" : "anything unusual or noteworthy",
              "onsets_adjust" : ( 0, 0, 0, 0 ),
              "run_st_mot_order" : ( 1, 2, 3, 0 )
            }

Pre-processing
--------------

Most of the pre-processing is done with the command ``glass_coherence_preproc``.
For help on using this script, run::

    glass_coherence_preproc --help

Typical usage is::

    glass_coherence_preproc sXXXX stage

where ``sXXXX`` is the subject ID and ``stage`` is the preprocessing stage (see below).

The stages are as follows:

Conversion
~~~~~~~~~~

Converts from the raw scanner format to a set of 4D NIFTI files::

    glass_coherence_preproc sXXXX convert

Correction
~~~~~~~~~~

Applies a motion and slice-timing correction procedure::

    glass_coherence_preproc sXXXX correct

*N.B. This stage takes quite a while...*

Fieldmaps
~~~~~~~~~

Prepares the fieldmaps::

    glass_coherence_preproc SXXXX fieldmaps

Unwarping
~~~~~~~~~

Before running, need to have made a symbolic link in each functional run directory to that run's fieldmap. For example::

    ln -s ../../fmap/f1/sXXXX_glass_coherence_fmap_1-fmap.nii sXXXX_glass_coherence_run_1-fmap.nii

Then, to use the fieldmaps to unwarp the functional images to remove the spatial distortion::

    glass_coherence_preproc sXXXX undistort

To verify that the unwarping has worked correctly:

* Run ``fslview``.
* Load the original or corrected image from a given run.
* Add the magnitude image from the fieldmap as an overlay.
* Notice the geometric distortions in the functional data.
* Add the undistorted image as an overlay, and hide the uncorrected image.
* Toggle the visibility of the undistorted image, and verify that the geometry now aligns well with that of the fieldmap's magnitude image.

ROI to images
~~~~~~~~~~~~~

Converts the raw ROI files from mrLoadRet into NIFTI masks::

    glass_coherence_preproc SXXXX roi-img

Coregistration
~~~~~~~~~~~~~~

The anatomical and ROI images are in a completely different space to the functionals, so they need to be coregistered.

The automatic FSL tools are *horrible* at doing this coregistration (in my experience), so we need to do it more manually using SPM.

Rough alignment
^^^^^^^^^^^^^^^

The coregistration algorithm is helped enormously if the images are in rough world-space alignment before it begins.

#. In SPM, click ``Display`` and select the mean functional image.
#. Place the crosshairs over a prominent landmark, such as the furthest posterior region of the occipital lobes. Note down the 3 values in the ``mm`` box.
#. Click ``Display`` again, this time selecting the anatomical image.
#. Place the crosshairs over the same landmark as was used in the functionals, and again note the 3 values in the ``mm`` box.
#. Subtract (element-wise) the anatomical ``mm`` values from the functional ``mm`` values, and use the output to populate the ``right``, ``forward``, and ``up`` fields.
#. To check your calculations, change the ``mm`` field to match what it was for the functional and the crosshairs should move to the same landmark.
#. Click ''Reorient images'' and select the anatomical **and the ROI images**.

Coregistration
^^^^^^^^^^^^^^

#. In SPM, click ``Coregister (Estimate & Reslice)``.
#. As the ``Reference image``, select the mean functional image.
#. As the ``Images to reslice``, select the anatomical image.
#. As the ``Other images``, select all the ROI images.
#. Under ``Reslice options``, change ``Interpolation`` to ``Nearest neighbour`` and ``Filename prefix`` to ``rs``.
#. Under ``File``, click ``Save batch`` and call it ``coreg.mat`` under the ``anat`` directory.
#. Click on the play icon to set it running.

Verification
^^^^^^^^^^^^

To check that the coregistration has performed well:

#. In SPM, click ``Check reg``.
#. Select the mean functional image first, and then the (unresliced) anatomical image.
#. Click around some prominent landmarks and check that the two images are in register.

ROI preparation
~~~~~~~~~~~~~~~

Converts the ROI image masks to a set of coordinates::

    glass_coherence_preproc sXXXX roi

Voxel timecourse extraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Extracts voxel timecourses for each voxel in each ROI::

    glass_coherence_preproc sXXXX vtc

Voxel culling
~~~~~~~~~~~~~

Removes voxels that have high mean-normalised variance::

    glass_coherence_preproc sXXXX vtc-cull

Timecourse averaging
~~~~~~~~~~~~~~~~~~~~

Averages over the voxels in each ROI::

    glass_coherence_preproc sXXXX vtc-avg

Design matrix
~~~~~~~~~~~~~

Forms the session design matrix from the run sequence information::

    glass_coherence_preproc sXXXX design


Subject-level analysis
----------------------

Most of the subject-level analysis is done with the command ``glass_coherence_subj_analysis``
For help on using this script, run::

    glass_coherence_subj_analysis --help

Typical usage is::

    glass_coherence_subj_analysis sXXXX stage

where ``sXXXX`` is the subject ID and ``stage`` is the preprocessing stage (see below).

The stages are as follows:

GLM
~~~

Fits a GLM to each ROI's timecourse::

    glass_coherence_subj_analysis sXXXX glm

Amplitude
~~~~~~~~~

This converts the GLM output to a single estimate per condition for a given subject::

    glass_coherence_subj_analysis sXXXX amp


Analysis datafiles
==================

The pre-processing / analysis pipeline produces the following files:

design
  ( 248 volumes, 4 conditions, 4 runs ) boolean array.
  This shows the onset of events for each condition.
  The matrix has been adjusted, where necessary, for delayed triggers.

task_info
  ( time, [ target, response, stim condition ], 4 runs ) integer array

coords-ROI
  ( 3 axes, n voxels ) array of coordinate locations.

coords_sel-ROI
  ( 3 axes, n voxels ) array of coordinate locations, *after* voxel selection.

vtc-ROI
  ( 248 volumes, 4 runs, n voxels ) array of BOLD signals.

vtc_sel-ROI
  ( 248 volumes, 4 runs, n voxels ) array of BOLD signals.
  This is as above but *after* voxel selection.

vtc_avg-ROI
  ( 248 volumes, 4 runs ) array of BOLD signals.

glm_beta-ROI
  ( n hrf, 4 conds, 4 runs ) array of beta values, in percent signal change units.

amp-ROI
  ( 4 conds ) vector of estimated HRF amplitudes, in percent signal change units.
