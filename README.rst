.. highlight:: tcsh

====================================================
Analysis code for natural scenes aperture experiment
====================================================

Requirements
============

- Python >= 2.7
- numpy
- matplotlib
- scipy
- fmri_tools (`http://www.bitbucket.org/djmannion/fmri_tools <http://www.bitbucket.org/djmannion/fmri_tools/>`_)


Processing stages
=================

Prepare the filesystem
----------------------

1. Make the subject's directory structure::

    mkdir -p sXXXX/{analysis/{exp,loc},fmap/f01,func/run{01,02,03,04,05,06,07,08,09,10,11,12},logs,rois,reg}

2. Copy the subject's runtime logfiles to the ``logs`` directory.

3. Make symlinks named ``raw`` in each functional run directory (exp and loc) that link to the location of its associated raw DICOM directory::

    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/MR-such_and_such raw

4. Similarly, make symlinks named ``mag-raw`` and ``ph-raw`` in each fieldmap directory that link to the locations of the fieldmap acquisition::

    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/MR-SEyada mag-raw
    ln -s /labs/olmanlab/DICOM/YYYYMMDD/sXXXX/PH-SEyada ph-raw

5. Make a local copy of the AFNI/SUMA base anatomicals (original and skull-stripped) that we can use for alignment::

    3dcopy \
       {$SUBJECTS_DIR}/{$SUBJ_ID}/SUMA/{$SUBJ_ID}_SurfVol+origK \
       reg/{$SUBJ_ID}_anat.nii


Update the experiment information file
--------------------------------------

Edit ``get_subj_conf`` within ``ns_aperture/config.py`` and add the new subject's information.

For example:

.. code-block:: python

    sXXXX = { "subj_id" : "sXXXX",
              "acq_date" : "YYYYMMDD",
              "n_runs" : 10,
              "n_loc_runs" : 2,
              "n_fmaps" : 1,
              "comments" : "anything unusual or noteworthy",
              "run_st_mot_order" : ( ( 7, "exp" ),
                                     ( 8, "exp" ),
                                     ( 9, "exp" ),
                                     ( 10, "exp" ),
                                     ( 1, "loc" ),
                                     ( 2, "loc" ),
                                     ( 1, "exp" ),
                                     ( 2, "exp" ),
                                     ( 3, "exp" ),
                                     ( 4, "exp" ),
                                     ( 5, "exp" ),
                                     ( 6, "exp" )
                                   ),
             "node_k" : { "lh" : 100000,
                          "rh" : 110000
                        }
            }

.. note::
   ``node_k`` can be found by viewing the subject's pial ASCII surface and noting the first number on the second row.


Pre-processing
--------------

Most of the pre-processing is done with the command ``ns_aperture_proc``.
For help on using this script, run::

    ns_aperture_proc --help

Typical usage is::

    ns_aperture_proc sXXXX stage

where ``sXXXX`` is the subject ID and ``stage`` is the preprocessing stage (see below).

The stages are as follows:

Conversion
~~~~~~~~~~

Converts from the raw scanner format to a set of 4D NIFTI files::

    ns_aperture_proc sXXXX convert

After execution, open up each NIFTI file and inspect for image quality and look at the summary image to see how much movement there was.


Correction
~~~~~~~~~~

Applies a slice-timing and motion correction procedure::

    ns_aperture_proc sXXXX correct

After execution, open up the summary NIFTI file to check that most of the motion has been removed.
You can also inspect the saved motion correction estimates to see how much movement there was.


Fieldmaps
~~~~~~~~~

Prepares the fieldmap::

    ns_aperture_proc SXXXX fieldmap


Unwarping
~~~~~~~~~

Before running, need to make a symbolic link in each functional run directory to that run's fieldmap. For example::

    ln -s ../../fmap/f01/sXXXX_ns_aperture_fmap_01-fmap.nii sXXXX_ns_aperture_run_01-fmap.nii

Then, to use the fieldmaps to unwarp the functional images to remove the spatial distortion::

    ns_aperture_proc sXXXX undistort

To verify that the unwarping has worked correctly:

* Run ``fslview``.
* Load the original or corrected image from a given run.
* Add the magnitude image from the fieldmap as an overlay.
* Notice the geometric distortions in the functional data.
* Add the undistorted image as an overlay, and hide the uncorrected image.
* Toggle the visibility of the undistorted image, and verify that the geometry now aligns well with that of the fieldmap's magnitude image.

Also, look at the session summary image produced and make sure that all looks good across the session.


Trim
~~~~

Removes timepoints from the start and/or end of each timeseries, as specified in the config::

    ns_aperture_proc SXXXX trim


Coregistration
~~~~~~~~~~~~~~

Follow the procedure described `here <http://visual-localiser-analysis-notes.readthedocs.org/en/latest/func.html#coregister-base-anatomy-to-functional-session>`__, substituting for the ``surf_reg`` command::

    ns_aperture_proc sXXXX surf_reg


Volume to surface
~~~~~~~~~~~~~~~~~

Projects the functional images to the cortical surface::

    ns_aperture_proc sXXXX vol_to_surf


Design preparation
~~~~~~~~~~~~~~~~~~

Computes the experimental design from the logfiles::

    ns_aperture_proc sXXXX design_prep


Subject-level analysis
----------------------

Localiser analysis
~~~~~~~~~~~~~~~~~~

Runs a GLM on the localiser data, extracts ``q`` (FDR) values, and creates a thresholded ROI mask::

    ns_aperture_proc sXXXX loc_glm


Experiment analysis
~~~~~~~~~~~~~~~~~~~

Runs a GLM on the experiment data::

    ns_aperture_proc sXXXX exp_glm


Datafile list
=============

Pre-processing
--------------


Subject-level analysis
----------------------


Group-level analysis
--------------------
