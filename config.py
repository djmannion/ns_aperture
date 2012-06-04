"""Configuration for the natural scenes through apertures fMRI experiment.
"""

from __future__ import division

import os

import numpy as np

import fmri_tools.utils, fmri_tools.paths


def get_conf():
	"""Overall experiment configuration.

	Returns
	-------
	conf : dict, with items:
		exp : holds overall experiment configurations
		stim : stimulus configuration
		task : behavioural task configuration
		acq : acquisition configuration
		ana : analysis configuration
		preproc : preprocessing configuration

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
	scale_mode : string, { "none", "cd", "norm", "mean" }
		How to scale the images.
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
	loc_rev_rate_hz : float
		Reversal rate of the localiser stimulus, in hertz.

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


def _get_exp_conf( tr_s = 2.0 ):
	"""Gets the experiment configuration.

	Specifies the experiment configuration and parameters for the natural scenes
	aperture fMRI experiment.

	This shouldn't be called directly.

	Parameters
	----------
	tr_s : float
		TR time, in seconds.

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
	rej_start_blocks : scalar int
		How many blocks to discard at the start of a run.
	rej_end_blocks : scalar int
		How many blocks to discard at the end of a run.
	hrf_corr_vol : scalar int
		Number of volumes to compensate for the HRF delay.
	n_vols_per_blk : scalar int
		Number of volumes per block.
	n_vols_per_run : scalar int
		Number of volumes per run.
	run_range_st : scalar int
		Starting index for desired run volume range.
	run_range_end : scalar int
		Ending index for desired run volume range.
	run_range : numpy array of integers
		Sequence of volume indices corresponding to valid run volumes.
	run_range_hrf_corr : numpy array of integers
		Sequence of volume indices corresponding to valid run volumes after HRF
		compensation.
	n_valid_vols_per_run : scalar int
		Number of volumes per run after trimming at the start and end.

	Notes
	-----
	* Returns are contained within a dictionary.

	"""

	exp_conf = {}

	exp_conf[ "id" ] = "ns_aperture"

	exp_conf[ "n_runs" ] = 10
	exp_conf[ "n_blocks" ] = 18
	exp_conf[ "block_len_s" ] = 16.0
	exp_conf[ "run_len_s" ] = exp_conf[ "n_blocks" ] * exp_conf[ "block_len_s" ]

	exp_conf[ "loc_n_runs" ] = 2
	exp_conf[ "loc_n_blocks" ] = 19
	exp_conf[ "loc_run_len_s" ] = ( exp_conf[ "loc_n_blocks" ] *
	                                exp_conf[ "block_len_s" ]
	                              )
	exp_conf[ "loc_pre_len_s" ] = 6

	exp_conf[ "loc_run_full_len_s" ] = ( exp_conf[ "loc_run_len_s" ] +
	                                     exp_conf[ "loc_pre_len_s" ]
	                                   )

	exp_conf[ "n_evt_per_block" ] = 12
	exp_conf[ "n_evt_per_run" ] = ( exp_conf[ "n_evt_per_block" ] *
	                                exp_conf[ "n_blocks" ]
	                              )
	exp_conf[ "evt_len_s" ] = ( exp_conf[ "block_len_s" ] /
	                            exp_conf[ "n_evt_per_block" ]
	                          )
	exp_conf[ "evt_stim_s" ] = 1.0

	exp_conf[ "rej_start_blocks" ] = 1
	exp_conf[ "rej_end_blocks" ] = 1

	exp_conf[ "n_vols_per_blk" ] = int( exp_conf[ "block_len_s" ] / tr_s )

	exp_conf[ "n_vols_per_run" ] = int( exp_conf[ "n_blocks" ] *
	                                    exp_conf[ "n_vols_per_blk" ]
	                                  )

	exp_conf[ "run_range_st" ] = int( exp_conf[ "rej_start_blocks" ] *
	                                  exp_conf[ "n_vols_per_blk" ]
	                                )
	exp_conf[ "run_range_end" ] = int( exp_conf[ "n_vols_per_run" ] -
	                                   exp_conf[ "rej_end_blocks" ] *
	                                   exp_conf[ "n_vols_per_blk" ]
	                                 )

	exp_conf[ "run_range" ] = np.arange( exp_conf[ "run_range_st" ],
	                                     exp_conf[ "run_range_end" ]
	                                   )

	exp_conf[ "n_valid_vols_per_run" ] = len( exp_conf[ "run_range" ] )

	exp_conf[ "n_valid_blocks" ] = ( exp_conf[ "n_blocks" ] -
	                                 exp_conf[ "rej_start_blocks" ] -
	                                 exp_conf[ "rej_end_blocks" ]
	                               )

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
	slice_axis : int
		the axis in the functional data that represents the inplanes.
	slice_acq_dir : { -1, +1 }
		whether the slices were acquired in the same order as represented in
		the array (1) or in descending order (-1).
	phase_encode_dir : { "x+", "y+", "z+", "x-", "y-", "z-" }
		FSLs unwarping code needs to know the phase encode direction (x|y|z) and
		polarity (-|+).

	Notes
	-----
	* Returns are contained within a dictionary.

	"""

	acq_conf = {}

	acq_conf = { "monitor_name" : "UMN_7T",
	             "tr_s" : 2.0,
	             "delta_te_ms" : 1.02,
	             "dwell_ms" : 0.325,
	             "slice_order" : fmri_tools.utils.get_slice_order( 36 ),
	             "slice_axis" : 2,
	             "slice_acq_dir" : 1,
	             "ph_encode_dir" : "y-",
	           }

	return acq_conf


def _get_preproc_conf():
	"""Gets the preprocessing configuration.

	Specifies the configuration for pre-processing the Glass pattern coherence
	fMRI data.

	Returns
	-------

	Notes
	-----
	* Return values are contained within a dictionary

"""

	slice_axis = 2

	slice_acq_dir = 1

	st_correct = False

	phase_encode_dir = "y-"

	roi_ax = ( 2, 1, 0 )

	roi_ax_order = ( 1, -1, -1 )

	# assemble the dictionary
	preproc_conf = { "slice_axis" : slice_axis,
	                 "slice_acq_dir" : slice_acq_dir,
	                 "phase_encode_dir" : phase_encode_dir,
	                 "roi_ax" : roi_ax,
	                 "roi_ax_order" : roi_ax_order,
	                 "st_correct" : st_correct
	               }

	return preproc_conf


def _get_analysis_conf():
	"""Gets the parameters for the fMRI analysis.

	Returns
	-------
	rois : tuple of strings
		Name of each ROI to be analysed.
	roi_ax : tuple of ints
		mapping from ROI axes to image axes. For example, ( 2, 1, 0 ) maps ROI
		(x,y,z) to image (z,y,x).
	roi_ax_order : tuple of { -1, +1 }
		the order when we map from ROI to image coordinates. -1 flips the axis
		order, while +1 preserves it.
	loc_p_thresh : float
		Probability threshold for the localiser analysis.
	poly_ord : integer >= 1
		Maximum Legendre polynomial order to use as nuisance regressors.

	Notes
	-----
	* Return values are within a dictionary.

	"""

	ana_conf = { "rois" : ( "V1",
	                        "V2",
	                        "V3",
	                        "V3AB",
	                        "hV4"
	                      ),
	             "roi_ax" : ( 2, 1, 0 ),
	             "roi_ax_order" : ( 1, -1, -1 ),
	             "loc_p_thresh" : 0.01,
	             "poly_ord" : 4,
	           }

	return ana_conf


def get_subj_conf( subj_id = None ):
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

	s1000 = { "subj_id" : "s1000",
	          "acq_date" : "20120601",
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

	subj_conf = { "s1000" : s1000,
	            }

	if subj_id is None:
		return subj_conf
	else:
		return subj_conf[ subj_id ]


def get_study_paths():
	"""Get the path structure for the study"""

	base_dir = "/labs/olmanlab/Data7T/NatSceneAperture/"

	study_paths = fmri_tools.paths.get_study_paths( base_dir )

	return study_paths


def get_subj_paths( subj_id ):
	"""Get the path structure for a given subject"""

	study_paths = get_study_paths()

	subj_conf = get_subj_conf( subj_id )

	study_conf = get_conf()

	subj_dir = os.path.join( study_paths[ "subj_dir" ], subj_id )

	# FUNCTIONALS
	func_dir = os.path.join( subj_dir, "func" )

	#   - experiment functionals
	func_exp_dir = os.path.join( func_dir, "exp" )

	func_paths = fmri_tools.paths.get_func_paths( func_exp_dir,
	                                              subj_id,
	                                              subj_conf[ "n_runs" ],
	                                              study_conf[ "exp" ][ "id" ],
	                                              sep_corr = False,
	                                              inc_xform_dirs = True
	                                            )

	#   - localiser functionals
	loc_id = "%s_loc" % study_conf[ "exp" ][ "id" ]
	func_loc_dir = os.path.join( func_dir, "loc" )

	loc_paths = fmri_tools.paths.get_func_paths( func_loc_dir,
	                                             subj_id,
	                                             subj_conf[ "n_loc_runs" ],
	                                             loc_id,
	                                             sep_corr = False,
	                                             inc_xform_dirs = True
	                                           )

	# SUMMARIES
	summ_paths = fmri_tools.paths.get_func_summ_paths( func_dir,
	                                                   subj_id,
	                                                   study_conf[ "exp" ][ "id" ],
	                                                   sep_corr = False
	                                                 )

	# FIELDMAPS
	fmap_dir = os.path.join( subj_dir, "fmap" )

	fmap_paths = fmri_tools.paths.get_fmap_paths( fmap_dir,
	                                              subj_id,
	                                              study_conf[ "exp" ][ "id" ],
	                                              subj_conf[ "n_fmaps" ]
	                                            )

	# ANATOMICALS
	anat_dir = os.path.join( subj_dir, "anat" )

	anat_paths = fmri_tools.paths.get_anat_paths( anat_dir,
	                                              subj_id,
	                                              study_conf[ "exp" ][ "id" ]
	                                            )

	# ROIS
	roi_dir = os.path.join( subj_dir, "roi" )
	roi_paths = fmri_tools.paths.get_roi_paths( roi_dir,
	                                            study_conf[ "ana" ][ "rois" ]
	                                          )

	subj_paths = { "func" : func_paths,
	               "loc" : loc_paths,
	               "summ" : summ_paths,
	               "fmap" : fmap_paths,
	               "anat" : anat_paths,
	               "roi" : roi_paths
	             }

	return subj_paths
