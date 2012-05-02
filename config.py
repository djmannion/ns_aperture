"""Configuration for the natural scenes through apertures fMRI experiment.
"""

from __future__ import division

import os
import string

import numpy as np

import fmri_tools.utils


def get_conf():
	"""Overall experiment configuration.

	Returns
	-------
	conf : dict, with items:
		exp : holds overall experiment configurations
		stim : stimulus configuration
		task : behavioural task configuration
		acq : acquisition configuration

	"""

	conf = { "exp" : _get_exp_conf(),
	         "stim" : _get_stim_conf(),
	         "task" : _get_task_conf(),
	         "acq" : _get_acq_conf(),
	         "ana" : _get_analysis_conf(),
	         "preproc" : _get_preproc_conf()
	       }

	return conf


def _get_stim_conf():
	"""Get the stimulus configuration.

	Specifies the stimulus configuration and parameters for the natural scenes
	aperture fMRI experiment.

	This shouldn't be called directly.

	Returns
	-------------
	db_path : string
		Path to the image database (the image files)
	im_dim : 2-item list
		( n_rows, n_cols ) dimensions of each (full) image
	im_deg_pp : scalar float
		Degrees per pixel in the images as acquired (ie native)
	patch_diam_deg : scalar float
		Diameter, in degrees visual angle, of each extracted image patch
	patch_ecc_deg : scalar float
		Eccenticity, in degrees visual angle, where each extracted patch will be
		shown
	patch_rect : 2-item list of 2-item list of 2-item list
		Coordinates of each patch, in the original image dimensions. Arranged as
		( patch 1, patch 2 ) -> ( top left, bottom right ) -> ( i_row, i_col )
	img_ids : numpy vector of ints
		Identifiers (one-based) for the set of candidate images.
	fix_rad_deg : float
		Radius of the fixation circle, in degrees visual angle.
	fix_col_inact : three-item tuple of floats
		Colour of inner fixation circle when task is inactive.
	fix_col_act : three-item tuple of floats
		Colour of inner fixation circle when task is active.
	loc_sf_cpd : float
		Spatial frequency of the localiser stimulus, in cycles per degree

	Notes
	-----
	* Returns are contained within a dictionary.

	"""

	stim_conf = {}

	stim_conf[ "db_path" ] = "../im_db"

	stim_conf[ "im_dim" ] = ( 1024, 1536 )
	stim_conf[ "im_deg_pp" ] = 1.0 / 60.0

	stim_conf[ "patch_diam_deg" ] = 4.0
	stim_conf[ "patch_ecc_deg" ] = 3.0

	patch_diam_pix = np.round( stim_conf[ "patch_diam_deg" ] *
	                           ( 1.0 / stim_conf[ "im_deg_pp" ] )
	                         )

	patch_r_pix = patch_diam_pix / 2

	patch_ecc_pix = np.round( stim_conf[ "patch_ecc_deg" ] *
	                          ( 1.0 / stim_conf[ "im_deg_pp" ] )
	                        )

	patch_centre_pix = [ ( ( stim_conf[ "im_dim" ][ 0 ] / 2 ),
	                       ( stim_conf[ "im_dim" ][ 1 ] / 2 +
	                         offset * patch_ecc_pix
	                       )
	                     )
	                     for offset in ( -1, +1 )
	                   ]

	stim_conf[ "patch_rect" ] = [ ( ( c_pix[ 0 ] - patch_r_pix, # top
	                                  c_pix[ 1 ] - patch_r_pix  # left
	                                ),
	                                ( c_pix[ 0 ] + patch_r_pix, # bottom
	                                  c_pix[ 1 ] + patch_r_pix  # right
	                                )
	                              )
	                              for c_pix in patch_centre_pix
	                            ]

	stim_conf[ "scale_mode" ] = "mean"

	# hardcode the desired image ids
	stim_conf[ "img_ids" ] = np.array( (    6,   12,  504,   58,   72,   90,
	                                       98,  110,  181,  230,  285,  342,
	                                      400,  460,  474,  482,  485,  561,
	                                      567,  639,  664,  727,  730,  737,
	                                      752,  756,  821, 1068, 1083, 1087,
	                                     1117, 1166, 1234, 1270, 1398, 1402,
	                                     1435, 1448, 1466, 1467, 1471, 1472,
	                                     1489, 1501, 1520, 1523, 1526, 1531,
	                                     1533, 1537, 1543, 1544, 1559, 1565,
	                                     1566, 1568, 1570, 1574, 1575, 1578,
	                                     1579, 1580, 1583, 1584, 1601, 1611,
	                                     1614, 1618, 1626, 1629, 1632, 1639,
	                                     1647, 1654, 1656, 1657, 1663, 1676,
	                                     1679, 1933, 1950, 1957, 2066, 2260,
	                                     2329, 2339, 2360, 2429, 2745, 2889,
	                                     2994, 2997, 3024, 3035, 3057, 3061,
	                                     3065, 3086, 3099, 3194, 3261, 3266,
	                                     3276, 3287, 3297, 3302, 3325, 3336
	                                  )
	                                )


	stim_conf[ "fix_rad_deg" ] = 0.075
	stim_conf[ "fix_col_inact" ] = ( -1, -1, -1 )
	stim_conf[ "fix_col_act" ] = ( -1, 1, -1 )

	stim_conf[ "loc_sf_cpd" ] = 2.0
	stim_conf[ "loc_rev_rate_hz" ] = 2

	return stim_conf


def _get_exp_conf():
	"""Gets the experiment configuration.

	Specifies the experiment configuration and parameters for the natural scenes
	aperture fMRI experiment.

	This shouldn't be called directly.

	Returns
	-------
	n_runs : scalar integer
		Number of experiment runs per session.
	n_blocks : scalar integer
		Number of blocks per run.
	block_len_s : scalar float
		Length of each block, in seconds.
	run_dur_s : scalar float
		Length of each run, in seconds.
	n_evt_per_block : scalar integer
		Number of events per block.
	n_evt_per_run : scalar integer
		Number of events per run.
	evt_len_s : scalar float
		Length of each 'event' within a block, in seconds.
	evt_stim_s : scalar float
		Length of stimulus presentation within each event.

	Notes
	-----
	* Returns are contained within a dictionary.

	"""

	exp_conf = {}

	exp_conf[ "n_runs" ] = 10

	exp_conf[ "n_blocks" ] = 18

	exp_conf[ "block_len_s" ] = 16.0

	exp_conf[ "run_len_s" ] = exp_conf[ "n_blocks" ] * exp_conf[ "block_len_s" ]

	exp_conf[ "n_evt_per_block" ] = 12

	exp_conf[ "n_evt_per_run" ] = ( exp_conf[ "n_evt_per_block" ] *
	                                exp_conf[ "n_blocks" ]
	                              )

	exp_conf[ "evt_len_s" ] = exp_conf[ "block_len_s" ] / exp_conf[ "n_evt_per_block" ]

	exp_conf[ "evt_stim_s" ] = 1.0

	exp_conf[ "rej_start_blocks" ] = 1
	exp_conf[ "rej_end_blocks" ] = 1

	exp_conf[ "hrf_corr_vol" ] = 2

	return exp_conf


def _get_task_conf():
	"""Gets the task configuration.

	Specifies the configuration for the behavioural task in the natural scenes
	aperture fMRI experiment.

	This shouldn't be called directly.

	Returns
	-------
	p : scalar float
		Probability that a given event will be a task trial.

	Notes
	-----
	* Returns are contained within a dictionary.

	"""

	task_conf = {}

	task_conf[ "p" ] = 0.35

	task_conf[ "evt_on_s" ] = 0.8
	task_conf[ "evt_off_s" ] = 1.2

	return task_conf


def _get_acq_conf():
	"""Get the acquisition configuration.

	Specifies the acquisition configuration and parameters for the
	natural scenes aperture fMRI experiment.

	Returns
	-------
	monitor_name : string
		monitor configuration name.
	tr_s : float
		time-to-repetition (TR), in seconds.
	delta_te_ms : float
		echo time differences for the fieldmaps, in milliseconds.
	dwell_ms : float
		dwell time, in milliseconds.
	slice_order : array of int
		slice acqusition indices, where 0 is the first slice.

	Notes
	-----
	* Returns are contained within a dictionary.

	"""

	acq_conf = {}

	acq_conf[ "monitor_name" ] = "UMN_7T"

	acq_conf[ "tr_s" ] = 2.0

	acq_conf[ "delta_te_ms" ] = 1.02

	acq_conf[ "dwell_ms" ] = 0.325

	acq_conf[ "slice_order" ] = fmri_tools.utils.get_slice_order( 36 )

	return acq_conf

def _get_preproc_conf():
	"""Gets the preprocessing configuration.

	Specifies the configuration for pre-processing the Glass pattern coherence
	fMRI data.

	Returns
	-------
	slice_axis : int
		the axis in the functional data that represents the inplanes.
	slice_acq_dir : { -1, +1 }
		whether the slices were acquired in the same order as represented in
		the array (1) or in descending order (-1).
	phase_encode_dir : { "x+", "y+", "z+", "x-", "y-", "z-" }
		FSLs unwarping code needs to know the phase encode direction (x|y|z) and
		polarity (-|+).
	roi_ax : tuple of ints
		mapping from ROI axes to image axes. For example, ( 2, 1, 0 ) maps ROI
		(x,y,z) to image (z,y,x).
	roi_ax_order : tuple of { -1, +1 }
		the order when we map from ROI to image coordinates. -1 flips the axis
		order, while +1 preserves it.

	Notes
	-----
	* Return values are contained within a dictionary

"""

	slice_axis = 2

	slice_acq_dir = 1

	phase_encode_dir = "y-"

	roi_ax = ( 2, 1, 0 )

	roi_ax_order = ( 1, -1, -1 )

	# assemble the dictionary
	preproc_conf = { "slice_axis" : slice_axis,
	                 "slice_acq_dir" : slice_acq_dir,
	                 "phase_encode_dir" : phase_encode_dir,
	                 "roi_ax" : roi_ax,
	                 "roi_ax_order" : roi_ax_order
	               }

	return preproc_conf


def _get_analysis_conf():
	"""Gets the parameters for the fMRI analysis.

	Returns
	-------
	rois : tuple of strings
		Name of each ROI to be analysed.
	hrf_len_vol : int, > 0
		Length of the HRF to estimate, in volumes.
	cull_prop : float, [0,1]
		Proportion of voxels in a given ROI that are culled based on their sorted
		mean-normalised variance.
	poly_ord : integer >= 1
		Maximum Legendre polynomial order to use as nuisance regressors.
	task_n_bins : integer >= 1
		Number of time windows (bins) to evaluate task performance.
	task_samp_rate_s : float
		The width of each bin, in seconds, in evaluating task performance.

	Notes
	-----
	* Return values are within a dictionary.

	"""

	rois = ( "V1",
	         "V2",
	         "V3",
	         "V3AB",
	         "hV4"
	       )

	loc_p_thresh = 0.01

	poly_ord = 4

	analysis_conf = { "rois" : rois,
	                  "loc_p_thresh" : loc_p_thresh,
	                  "poly_ord" : poly_ord,
	                }

	return analysis_conf


def get_subj_conf():
	"""Gets the configuration info for each subject.

	Returns
	-------
	subj_id : string
		Subject ID, in the Olman lab system.
	acq_date : string
		Acquisition date, in YYYYMMDD format.
	comments : string
		Any comments about the scanning session.

	Notes
	-----
	* Return values are within a dictionary, which is itself within a dictionary
	  indexed by subject ID.

	"""

	s1021 = { "subj_id" : "s1021",
	          "acq_date" : "20120430",
	          "n_runs" : 10,
	          "n_loc_runs" : 2,
	          "n_fmaps" : 1,
	          "run_st_mot_order" : ( ( 7, "func" ),
	                                 ( 8, "func" ),
	                                 ( 9, "func" ),
	                                 ( 10, "func" ),
	                                 ( 1, "loc" ),
	                                 ( 2, "loc" ),
	                                 ( 1, "func" ),
	                                 ( 2, "func" ),
	                                 ( 3, "func" ),
	                                 ( 4, "func" ),
	                                 ( 5, "func" ),
	                                 ( 6, "func" )
	                               ),
	          "comments" : ""
	        }

	subj_conf = { "s1021" : s1021,
	            }

	return subj_conf


def get_exp_paths():
	"""
	"""

	# hard code
	BASE_PATH = "/home/dmannion/NatSceneAperture/"

	paths = {}

	paths[ "dir" ] = BASE_PATH

	paths[ "grp_dir" ] = os.path.join( paths[ "dir" ],
	                                   "group_data"
	                                 )

	paths[ "grp_task" ] = os.path.join( paths[ "grp_dir" ],
	                                    "grp_task.mat"
	                                  )

	return paths


def get_subj_paths( subj_id ):
	"""Gets the filesystem path and file structure for the experiment data

	Parameters
	----------
	subj_id : string
		Subject identification number, eg. s1021

	Returns
	-------
	A dictionary with the following path info
	[ "exp" ]
		Paths at the experiment level
	[ "subj" ]
		Paths at the subject level
	[ "dir" ]
		Base directory for the subject.
	[ "func" ]
		Paths for the functional data.
	[ "fmap" ]
		Paths for the fieldmap data.
	[ "anat" ]
		Paths for the anatomical data.
	[ "roi" ]
		Paths for the ROI data.
	[ "analysis" ]
		Paths for the data analysis.
	[ "task" ]
		Paths for the task data.
	[ "design" ]
		Path to the design matrix.
	"""

	def get_task_paths( subj_conf, subj_dir ):
		""" Returns the paths for the fixation task

		Parameters
		----------
		subj_conf: dict
			Subject configuration directory
		subj_dir : string
			Path to subject base directory

		Returns
		-------
		task : dict of strings
			Dictionary with the following path info:
			[ "dir" ] :
				Base directory for the task.
			[ "task_info" ] :
				Base path for the task information data.
			[ "task_resp" ] :
				Base path for the task response data.
		"""

		task = {}

		task[ "dir" ] = os.path.join( subj_dir, "log" )

		task[ "task_info" ] = os.path.join( task[ "dir" ],
		                                    "%s_ns_aperture_info.npy" % (
		                                    subj_conf[ "subj_id" ] )
		                                  )

		return task


	def get_analysis_paths( subj_dir ):
		""" Returns the paths for the analysis files

		Parameters
		----------
		subj_conf: dict
			Subject configuration directory
		subj_dir : string
			Path to subject base directory

		Returns
		-------
		analysis : dict of strings
			Dictionary with the following path info:
			[ "dir" ] :
				Base directory for the analysis.
			[ "vtc" ] :
				Base path for the voxel timecourses of each ROI.
			[ "coords" ] :
				Base path for the coordinates of each ROI.
		"""

		analysis = {}

		analysis[ "dir" ] = os.path.join( subj_dir, "analysis" )

		analysis[ "vtc" ] = os.path.join( analysis[ "dir" ], "vtc" )
		analysis[ "loc_vtc" ] = os.path.join( analysis[ "dir" ], "loc_vtc" )

		analysis[ "loc_stat" ] = os.path.join( analysis[ "dir" ], "loc_stat" )

		analysis[ "vtc_sel" ] = os.path.join( analysis[ "dir" ], "vtc_sel" )
		analysis[ "vtc_avg" ] = os.path.join( analysis[ "dir" ], "vtc_avg" )

		analysis[ "coords" ] = os.path.join( analysis[ "dir" ], "coords" )
		analysis[ "coords_sel" ] = os.path.join( analysis[ "dir" ], "coords_sel" )

		analysis[ "glm_beta" ] = os.path.join( analysis[ "dir" ], "glm_beta" )
		analysis[ "glm_tc" ] = os.path.join( analysis[ "dir" ], "glm_tc" )
		analysis[ "glm_gof" ] = os.path.join( analysis[ "dir" ], "glm_gof" )

		analysis[ "amp" ] = os.path.join( analysis[ "dir" ], "amp" )

		return analysis


	def get_func_paths( subj_conf, subj_dir ):
		""" Returns the paths for the functional acquisitions

		Parameters
		----------
		subj_conf : dict
			Subject configuration dictionary
		subj_dir : string
			Path to subject base directory

		Returns
		-------
		func : dict of strings
			Dictionary with the following path info:
			[ "dir " ] :
				Base directory for the functional data.
			[ "dirs" ] :
				Directory of each functional run.
			[ "dcm_dirs" ] :
				Directory of each functional run's raw DICOM directory.
			[ "orig" ], [ "corr" ], [ "uw" ]
				Filenames for each functional run's original, corrected, and unwarped
				NII files.
			[ "fmap" ] :
				Filename of each functional run's fieldmap.
			[ "mean" ] :
				Filename of an image containing the mean of all functional data.
			[ "motion_estimates" ] :
				Filename of an array containing the estimated motion parameters.

		"""

		func = {}

		# base functional directory
		func[ "dir" ] = os.path.join( subj_dir, "func" )

		# directory of each functional run
		func[ "dirs" ] = [ os.path.join( func[ "dir" ],
		                                 "run%02d" % ( i_run + 1 )
		                               )
		                   for i_run in xrange( subj_conf[ "n_runs" ] )
		                 ]

		# directory of each functional run's raw dicoms
		func[ "dcm_dirs" ] = [ os.path.join( func_dir,
		                                     "raw"
		                                   )
		                       for func_dir in func[ "dirs" ]
		                     ]

		# path to each functional run's 4D nifti file
		func[ "orig" ] = [ os.path.join( func_dir,
		                                 "%s_ns_aperture_run_%02d-orig" %
		                                 ( subj_conf[ "subj_id" ],
		                                   i_run + 1
		                                 )
		                               )
		                   for ( i_run, func_dir ) in enumerate( func[ "dirs" ] )
		                 ]

		# path to each functional run's 4D nifti file, slice-time and motion
		# corrected
		func[ "corr" ] = [ ( "%s-corr" % orig )
		                   for orig in func[ "orig" ]
		                 ]

		# ... and unwarped (distortion corrected)
		func[ "uw" ] = [ ( "%s-uw" % corr )
		                 for corr in func[ "corr" ]
		               ]

		# path to the fieldmap for each functional run
		func[ "fmap" ] = [ string.replace( orig, "-orig", "-fmap" )
		                   for orig in func[ "orig" ]
		                 ]

		# path to a mean image of all functional images
		func[ "mean" ] = os.path.join( func[ "dir" ],
		                               "%s_ns_aperture-mean" %
		                               subj_conf[ "subj_id" ]
		                             )

		# path to a mean image of all functional images
		func[ "summ" ] = os.path.join( func[ "dir" ],
		                               "%s_ns_aperture-summ" %
		                               subj_conf[ "subj_id" ]
		                             )

		# path to numpy file containing the motion estimates from the correction
		# procedure
		func[ "motion_estimates" ] = os.path.join( func[ "dir" ],
		                                           "motion_estimates.npy"
		                                         )

		return func

	def get_loc_paths( subj_conf, subj_dir ):
		""" Returns the paths for the functional acquisitions

		Parameters
		----------
		subj_conf : dict
			Subject configuration dictionary
		subj_dir : string
			Path to subject base directory

		Returns
		-------
		func : dict of strings
			Dictionary with the following path info:
			[ "dir " ] :
				Base directory for the functional data.
			[ "dirs" ] :
				Directory of each functional run.
			[ "dcm_dirs" ] :
				Directory of each functional run's raw DICOM directory.
			[ "orig" ], [ "corr" ], [ "uw" ]
				Filenames for each functional run's original, corrected, and unwarped
				NII files.
			[ "fmap" ] :
				Filename of each functional run's fieldmap.
			[ "mean" ] :
				Filename of an image containing the mean of all functional data.
			[ "motion_estimates" ] :
				Filename of an array containing the estimated motion parameters.

		"""

		loc = {}

		# base functional directory
		loc[ "dir" ] = os.path.join( subj_dir, "loc" )

		# directory of each functional run
		loc[ "dirs" ] = [ os.path.join( loc[ "dir" ],
		                                "run%02d" % ( i_run + 1 )
		                              )
		                  for i_run in xrange( subj_conf[ "n_loc_runs" ] )
		                ]

		# directory of each functional run's raw dicoms
		loc[ "dcm_dirs" ] = [ os.path.join( loc_dir,
		                                    "raw"
		                                  )
		                      for loc_dir in loc[ "dirs" ]
		                    ]

		# path to each functional run's 4D nifti file
		loc[ "orig" ] = [ os.path.join( loc_dir,
		                                "%s_ns_aperture_loc_run_%02d-orig" %
		                                 ( subj_conf[ "subj_id" ],
		                                   i_run + 1
		                                 )
		                               )
		                  for ( i_run, loc_dir ) in enumerate( loc[ "dirs" ] )
		                ]

		# path to each functional run's 4D nifti file, slice-time and motion
		# corrected
		loc[ "corr" ] = [ ( "%s-corr" % orig )
		                   for orig in loc[ "orig" ]
		                 ]

		# ... and unwarped (distortion corrected)
		loc[ "uw" ] = [ ( "%s-uw" % corr )
		                for corr in loc[ "corr" ]
		              ]

		# path to the fieldmap for each functional run
		loc[ "fmap" ] = [ string.replace( orig, "-orig", "-fmap" )
		                  for orig in loc[ "orig" ]
		                ]

		return loc


	def get_fmap_paths( subj_conf, subj_dir ):
		""" Returns the paths for the fieldmap acquisitions

		Parameters
		----------
		subj_conf : dict
			Subject configuration dictionary
		subj_dir : string
			Path to subject base directory

		Returns
		-------
		fmap : dict of paths
			[ "dir " ] :
				Base directory for the fieldmap data.
			[ "dirs" ] :
				Directory of each fieldmap acquisition.
			[ "dcm_mag_dirs" ], [ "dcm_ph_dirs" ]:
				Directory of each fieldmap's magnitude and phase raw DICOM directories.
			[ "mag" ], [ "ph" ] :
				Filenames for each fieldmap's magnitude and phase NII files.
			[ "fmap" ] :
				Filename of each fieldmap's final fieldmap.

		"""

		fmap = {}

		# base directory for the fieldmaps
		fmap[ "dir" ] = os.path.join( subj_dir, "fmap" )

		# directories of each of the fieldmaps acquired in the session
		fmap[ "dirs" ] = [ os.path.join( fmap[ "dir" ],
		                                 "f%d" % ( i_fmap + 1 )
		                               )
		                   for i_fmap in xrange( subj_conf[ "n_fmaps" ] )
		                 ]

		# directories holding the raw magnitude images for each fieldmap
		fmap[ "dcm_mag_dirs" ] = [ os.path.join( fmap_dir,
		                                         "mag-raw"
		                                       )
		                           for fmap_dir in fmap[ "dirs" ]
		                         ]

		# ... raw phase image directories
		fmap[ "dcm_ph_dirs" ] = [ os.path.join( fmap_dir,
		                                        "ph-raw"
		                                      )
		                          for fmap_dir in fmap[ "dirs" ]
		                        ]

		# paths to the nifti file of each fieldmap's magnitude image
		fmap[ "mag" ] = [ os.path.join( fmap_dir,
		                                "%s_ns_aperture_fmap_%d-mag" %
		                                ( subj_conf[ "subj_id" ], i_fmap + 1 )
		                              )
		                  for ( i_fmap, fmap_dir ) in enumerate( fmap[ "dirs" ] )
		                ]

		# ... phase image
		fmap[ "ph" ] = [ os.path.join( fmap_dir,
		                               "%s_ns_aperture_fmap_%d-ph" %
		                               ( subj_conf[ "subj_id" ], i_fmap + 1 )
		                             )
		                 for ( i_fmap, fmap_dir ) in enumerate( fmap[ "dirs" ] )
		               ]

		# paths to the nifti file of each fieldmap
		fmap[ "fmap" ] = [ os.path.join( fmap_dir,
		                                 "%s_ns_aperture_fmap_%d-fmap" %
		                                 ( subj_conf[ "subj_id" ], i_fmap + 1 )
		                               )
		                   for ( i_fmap, fmap_dir ) in enumerate( fmap[ "dirs" ] )
		                 ]

		return fmap


	def get_design_paths( subj_conf, subj_dir ):
		""" Returns the paths for the session design.
		Parameters
		----------
		subj_conf : dict
			Subject configuration dictionary
		subj_dir : string
			Path to subject base directory

		Returns
		-------
		design : dict of paths
			[ "matrix" ] :
				Path to the design matrix.
			[ "run_seq" ] :
				Paths to the sequence info for each run.

		"""

		design = {}

		design[ "dir" ] = os.path.join( subj_dir, "log" )

		design[ "matrix" ] = os.path.join( subj_dir,
		                                   "design_matrix.npy"
		                                 )

		design[ "log" ] = os.path.join( design[ "dir" ],
		                                "%s_ns_aperture.mat" %
		                                subj_conf[ "subj_id" ]
		                              )

		return design


	def get_anat_paths( subj_conf, subj_dir ):
		""" Returns the paths for the anatomy

		Parameters
		----------
		subj_conf : dict
			Subject configuration dictionary
		subj_dir : string
			Path to subject base directory

		Returns
		-------
		anat : dict of paths
			[ "dir" ] :
				Base directory for the anatomical data.
			[ "anat" ] :
				Filename of the anatomical image.

		"""

		anat = {}

		# base anatomical directory
		anat[ "dir" ] = os.path.join( subj_dir, "anat" )

		# path to the anatomical nifti
		anat[ "anat" ] = os.path.join( anat[ "dir" ],
		                               "%s_anat" % subj_conf[ "subj_id" ]
		                             )

		return anat


	def get_roi_paths( subj_dir, roi_names ):
		""" Returns the paths for the ROIs

		Parameters
		----------
		subj_dir : string
			Path to subject base directory
		roi_names : list of strings
			Names of the ROIs

		Returns
		-------
		roi : dict of paths
			[ "dir" ] :
				Base directory for the ROI data.
			[ "mat" ] :
				Filenames for the ROI MAT files.
			[ "orig" ] :
				Filenames of the original ROI images.
			[ "rs" ] :
				Filenames of the resampled ROI images.

		"""

		roi = {}

		# base ROI directory
		roi[ "dir" ] = os.path.join( subj_dir, "roi" )

		# path to mat file for each ROI
		roi[ "mat" ] = [ os.path.join( roi[ "dir" ],
		                               "%s.mat" % roi_name
		                             )
		                 for roi_name in roi_names
		               ]

		# path to the nifti file of each ROI, in its original (anatomical) image
		# space
		roi[ "orig" ] = [ os.path.join( roi[ "dir" ],
		                               "%s_mask" % roi_name
		                              )
		                  for roi_name in roi_names
		                ]

		# path to the nifti file of each ROI in the same space as the functional
		# images
		roi[ "rs" ] = [ os.path.join( roi[ "dir" ],
		                              "rs_%s_mask" % roi_name
		                            )
		                for roi_name in roi_names
		              ]

		return roi


	analysis_conf = get_conf()[ "ana" ]

	subj_conf = get_subj_conf()[ subj_id ]

	paths = {}

	paths[ "exp" ] = get_exp_paths()

	paths[ "subj" ] = { "dir" : os.path.join( paths[ "exp" ][ "dir" ],
	                                          "subj_data",
	                                          subj_conf[ "subj_id" ]
	                                        )
	                  }

	paths[ "func" ] = get_func_paths( subj_conf, paths[ "subj" ][ "dir" ] )

	paths[ "loc" ] = get_loc_paths( subj_conf, paths[ "subj" ][ "dir" ] )

	paths[ "fmap" ] = get_fmap_paths( subj_conf, paths[ "subj" ][ "dir" ] )

	paths[ "anat" ] = get_anat_paths( subj_conf, paths[ "subj" ][ "dir" ] )

	paths[ "roi" ] = get_roi_paths( paths[ "subj" ][ "dir" ],
	                                analysis_conf[ "rois" ]
	                              )

	paths[ "design" ] = get_design_paths( subj_conf, paths[ "subj" ][ "dir" ] )

	paths[ "analysis" ] = get_analysis_paths( paths[ "subj" ][ "dir" ] )

	paths[ "task" ] = get_task_paths( subj_conf, paths[ "subj" ][ "dir" ] )

	return paths

