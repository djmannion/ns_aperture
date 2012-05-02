"""
Set of routines to pre-process the fMRI data for the Glass patterns coherence
fMRI experiment.
"""

from __future__ import division

import os.path as path
import itertools

import nipy
import numpy as np
import scipy.io

import fmri_tools.preproc, fmri_tools.utils

def convert( paths ):
	"""Converts the functionals and fieldmaps from dicom to nifti.

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subject_paths' in
		'glass_coherence.config'.

	"""

	# functionals
	func_names = [ path.split( func_path )[ 1 ]
	               for func_path in paths[ "func" ][ "orig" ]
	             ]

	map( fmri_tools.preproc.dcm_to_nii,
	     paths[ "func" ][ "dcm_dirs" ],  # dcm path
	     paths[ "func" ][ "dirs" ],  # nii path
	     func_names  # nii name
	   )

	# check that there hasn't been any mix-ups with the assignment of raw volumes
	# to runs
	# (need to add the extension)
	full_func_paths = [ "".join( [ func_path, ".nii" ] )
	                    for func_path in paths[ "func" ][ "orig" ]
	                  ]

	assert( fmri_tools.utils.files_are_unique( full_func_paths ) )

	# localisers
	loc_names = [ path.split( loc_path )[ 1 ]
	              for loc_path in paths[ "loc" ][ "orig" ]
	            ]

	map( fmri_tools.preproc.dcm_to_nii,
	     paths[ "loc" ][ "dcm_dirs" ],  # dcm path
	     paths[ "loc" ][ "dirs" ],  # nii path
	     loc_names  # nii name
	   )

	# check that there hasn't been any mix-ups with the assignment of raw volumes
	# to runs
	# (need to add the extension)
	full_loc_paths = [ "".join( [ loc_path, ".nii" ] )
	                   for loc_path in paths[ "loc" ][ "orig" ]
	                 ]

	assert( fmri_tools.utils.files_are_unique( full_loc_paths ) )

	# fieldmaps - magnitude
	mag_names = [ path.split( mag_path )[ 1 ]
	              for mag_path in paths[ "fmap" ][ "mag" ]
	            ]

	map( fmri_tools.preproc.dcm_to_nii,
	     paths[ "fmap" ][ "dcm_mag_dirs" ],  # dcm path
	     paths[ "fmap" ][ "dirs" ],  # nii path
	     mag_names  # nii name
	   )

	full_mag_paths = [ "".join( [ mag_path, ".nii" ] )
	                    for mag_path in paths[ "fmap" ][ "mag" ]
	                  ]

	assert( fmri_tools.utils.files_are_unique( full_mag_paths ) )

	# fieldmaps - phase
	ph_names = [ path.split( ph_path )[ 1 ]
	             for ph_path in paths[ "fmap" ][ "ph" ]
	           ]

	map( fmri_tools.preproc.dcm_to_nii,
	     paths[ "fmap" ][ "dcm_ph_dirs" ],  # dcm path
	     paths[ "fmap" ][ "dirs" ],  # nii path
	     ph_names  # nii name
	   )

	full_ph_paths = [ "".join( [ ph_path, ".nii" ] )
	                   for ph_path in paths[ "fmap" ][ "ph" ]
	                 ]

	assert( fmri_tools.utils.files_are_unique( full_ph_paths ) )


def st_motion_correct( paths, acq_conf, preproc_conf, subj_conf ):
	"""Performs slice-timing and motion correction

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'glass_coherence.config'.
	acq_conf : dict
		Experiment configuration, as returned by 'get_acq_conf' in
		'glass_coherence.config'.
	preproc_conf : dict
		Preprocessing configuration, as returned by 'get_preproc_conf' in
		'glass_coherence.config'.
	subj_conf : dict
		Subject configuration, as returned by 'get_subj_conf' in
		'glass_coherence.config', for this subject.

	"""

	# get the order of runs to pass to the correction algorithm
	# this is done because the algorithm realigns all to the first entry, which
	# normally corresponds to the first run acquired, but we might want them to
	# be aligned with a different run - one closer to a fieldmap, for example
	run_order = subj_conf[ "run_st_mot_order" ]

	# reorder the paths

	orig_paths = [ paths[ im_type ][ "orig" ][ i_run - 1 ]
	               for i_run, im_type in run_order
	             ]
	corr_paths = [ paths[ im_type ][ "corr" ][ i_run - 1 ]
	               for i_run, im_type in run_order
	             ]


	motion_est = fmri_tools.preproc.correct_st_motion( orig_paths,
	                                                   corr_paths,
	                                                   acq_conf[ "slice_order" ],
	                                                   acq_conf[ "tr_s" ],
	                                                   ( preproc_conf[ "slice_axis" ],
	                                                     preproc_conf[ "slice_acq_dir" ]
	                                                   )
	                                                 )

	np.save( paths[ "func" ][ "motion_estimates" ],
	         arr = motion_est
	       )


def fieldmaps( paths, acq_conf, subj_conf ):
	"""Prepare the fieldmaps.

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'glass_coherence.config'.
	acq_conf : dict
		Experiment configuration, as returned by 'get_acq_conf' in
		'glass_coherence.config'.
	subj_conf : dict
		Subject configuration, as returned by 'get_subj_conf' in
		'glass_coherence.config', for this subject.
	"""

	map( fmri_tools.preproc.make_fieldmap,
	     paths[ "fmap" ][ "mag" ],
	     paths[ "fmap" ][ "ph" ],
	     paths[ "fmap" ][ "fmap" ],
	     itertools.repeat( acq_conf[ "delta_te_ms" ],
	                       subj_conf[ "n_fmaps" ]
	                     )
	   )


def unwarp( paths, acq_conf, preproc_conf, subj_conf ):
	"""Uses the fieldmaps to unwarp the functional images and create a mean image
	of all the unwarped functional images.

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'glass_coherence.config'.
	acq_conf : dict
		Experiment configuration, as returned by 'get_acq_conf' in
		'glass_coherence.config'.
	preproc_conf : dict
		Preprocessing configuration, as returned by 'get_preproc_conf' in
		'glass_coherence.config'.
	subj_conf : dict
		Subject configuration, as returned by 'get_subj_conf' in
		'glass_coherence.config', for this subject.

	"""

	# perform the unwarping on the funcs
	map( fmri_tools.preproc.unwarp,
	     paths[ "func" ][ "corr" ],
	     paths[ "func" ][ "fmap" ],
	     paths[ "func" ][ "uw" ],
	     itertools.repeat( acq_conf[ "dwell_ms" ],
	                       subj_conf[ "n_runs" ]
	                     ),
	     itertools.repeat( preproc_conf[ "phase_encode_dir" ],
	                       subj_conf[ "n_runs" ]
	                     )
	   )

	# create a mean image of the unwarped data
	fmri_tools.preproc.mean_image( paths[ "func" ][ "uw" ],
	                               paths[ "func" ][ "mean" ]
	                             )

	# perform the unwarping on the localisers
	map( fmri_tools.preproc.unwarp,
	     paths[ "loc" ][ "corr" ],
	     paths[ "loc" ][ "fmap" ],
	     paths[ "loc" ][ "uw" ],
	     itertools.repeat( acq_conf[ "dwell_ms" ],
	                       subj_conf[ "n_loc_runs" ]
	                     ),
	     itertools.repeat( preproc_conf[ "phase_encode_dir" ],
	                       subj_conf[ "n_loc_runs" ]
	                     )
	   )

	# produce a summary image
	fmri_tools.preproc.gen_sess_summ_img( ( paths[ "func" ][ "uw" ] +
	                                        paths[ "loc" ][ "uw" ]
	                                      ),
	                                      paths[ "func" ][ "summ" ]
	                                    )


def make_roi_images( paths, ana_conf, preproc_conf ):
	"""Converts the ROI matlab files to nifti images in register with the
	subject's anatomical.

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'glass_coherence.config'.
	ana_conf : dict
		Analysis configuration, as returned by 'get_analysis_conf' in
		'glass_coherence.config'.
	preproc_conf : dict
		Preprocessing configuration, as returned by 'get_preproc_conf' in
		'glass_coherence.config', for this subject.

	"""

	n_rois = len( ana_conf[ "rois" ] )

	# python nifti loader needs the extension
	anat_path = "".join( [ paths[ "anat" ][ "anat" ],
	                       ".nii"
	                     ]
	                   )

	roi_paths = [ "".join( [ roi_path, ".nii" ] )
	              for roi_path in paths[ "roi" ][ "orig" ]
	            ]

	# convert the ROI mat files to nifti images
	map( fmri_tools.preproc.roi_to_nii,
	     paths[ "roi" ][ "mat" ],
	     itertools.repeat( anat_path, n_rois ),
	     roi_paths,
	     itertools.repeat( preproc_conf[ "roi_ax" ], n_rois ),
	     itertools.repeat( preproc_conf[ "roi_ax_order" ], n_rois )
	   )


def prepare_rois( paths, ana_conf ):
	"""Extracts and writes the coordinates of each ROI.

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'glass_coherence.config'.
	ana_conf : dict
		Analysis configuration, as returned by 'get_analysis_conf' in
		'glass_coherence.config'.

	"""

	for ( i_roi, roi_name ) in enumerate( ana_conf[ "rois" ] ):

		# load the ROI image
		roi = nipy.load_image( "".join( [ paths[ "roi" ][ "rs" ][ i_roi ],
		                                  ".nii"
		                                ]
		                              )
		                     ).get_data()

		# get rid of any NaNs
		roi[ np.isnan( roi) ] = 0

		# extract the coordinate list
		coords = np.nonzero( roi )

		# and save
		np.save( "".join( [ paths[ "analysis" ][ "coords" ],
		                    "-",
		                    roi_name
		                  ]
		                ),
		         arr = coords
		       )


def roi_vtc_cull( paths, ana_conf ):
	"""Culls voxels from each ROI based on their mean-normalised variance.

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'glass_coherence.config'.
	ana_conf : dict
		Analysis configuration, as returned by 'get_analysis_conf' in
		'glass_coherence.config'.

	"""

	for roi_name in ana_conf[ "rois" ]:

		# load the coordinates
		coords = np.load( "".join( [ paths[ "analysis" ][ "coords" ],
		                             "-", roi_name, ".npy"
		                           ]
		                         )
		                )

		# load the vtc
		vtc = np.load( "".join( [ paths[ "analysis" ][ "vtc" ],
		                          "-", roi_name, ".npy"
		                        ]
		                      )
		             )

		# do a quick check to make sure the coords have the expected shape
		assert( coords.shape[ 0 ] == 3 )
		assert( coords.shape[ 1 ] == vtc.shape[ -1 ] )

		n_voxels = coords.shape[ 1 ]

		# this is a logical matrix showing if the voxel passes the selection
		# criteria for each run
		retain = np.empty( ( vtc.shape[ 1:3 ] ) )

		# this is the number of voxels post-culling (for each run)
		cutoff_n = np.round( ( 1 - ana_conf[ "cull_prop" ] ) * n_voxels )

		for i_run in xrange( vtc.shape[ 1 ] ):

			# decision statistic is mean-normalised variance
			dec_stat = ( np.var( vtc[ :, i_run, : ], axis = 0 ) /
			             np.mean( vtc[ :, i_run, : ], axis = 0 )
			           )

			# sort in ascending order
			dec_stat_sorted = np.sort( dec_stat )

			# extract the critical stat value
			crit_stat = dec_stat_sorted[ cutoff_n ]

			retain[ i_run, : ] = dec_stat < crit_stat

		# to be retained, has to be retained for all runs
		vox_to_retain = np.all( retain, axis = 0 )

		# do the cullin'
		vtc = vtc[ :, :, vox_to_retain ]
		coords = coords[ :, vox_to_retain ]

		# and re-save
		np.save( "".join( [ paths[ "analysis" ][ "coords_sel" ],
		                           "-", roi_name, ".npy"
		                         ]
		                ),
		         coords
		       )

		# save the vtc
		np.save( "".join( [ paths[ "analysis" ][ "vtc_sel" ],
		                           "-", roi_name, ".npy"
		                         ]
		                ),
		         vtc
		       )

def avg_vtcs( paths, ana_conf ):
	"""Averages the timecourses over all the (selected) voxels in a ROI.

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'glass_coherence.config'.
	ana_conf : dict
		Analysis configuration, as returned by 'get_analysis_conf' in
		'glass_coherence.config'.

	"""

	for roi_name in ana_conf[ "rois" ]:

		# load the (post voxel selection) vtc
		vtc = np.load( "".join( [ paths[ "analysis" ][ "vtc_sel" ],
		                          "-", roi_name, ".npy"
		                        ]
		                      )
		             )

		# average over voxels
		vtc = np.mean( vtc, axis = 2 )

		# save
		np.save( "".join( [ paths[ "analysis" ][ "vtc_avg" ],
		                   "-", roi_name
		                  ]
		                ),
		         arr = vtc
		       )


def localiser_analysis( paths, exp_conf, ana_conf, acq_conf ):
	"""
	"""

	block_len_vol = int( exp_conf[ "block_len_s" ] / acq_conf[ "tr_s" ] )

	for roi_name in ana_conf[ "rois" ]:

		# load the (post voxel selection) vtc
		loc_vtc = np.load( "%s-%s.npy" % ( paths[ "analysis" ][ "loc_vtc" ],
		                                   roi_name
		                                 )
		                 )

		( n_vol, n_runs, n_voxels ) = loc_vtc.shape

		block_lut = np.empty( ( n_runs, int( n_vol / block_len_vol ), 2 ) )

		block_lut[ 0, :, 0 ] = np.arange( 0, n_vol, block_len_vol )
		block_lut[ 0, :, 1 ] = np.tile( [ 0, 1 ], int( n_vol / block_len_vol / 2 ) )

		block_lut[ 1, :, 0 ] = block_lut[ 0, :, 0 ]
		block_lut[ 1, :, 1 ] = np.tile( [ 1, 0 ], int( n_vol / block_len_vol / 2 ) )

		stat = np.empty( ( n_voxels, 2 ) )

		for i_voxel in xrange( n_voxels ):

			vox_data = [ [], [] ]

			for i_run in xrange( n_runs ):

				run_vtc = loc_vtc[ :, i_run, i_voxel ]
				run_vtc = fmri_tools.preproc.hp_filter( run_vtc, ana_conf[ "poly_ord" ] )[ 0 ]

				for i_blk in xrange( block_lut.shape[ 1 ] ):

					blk_range = np.arange( block_lut[ i_run, i_blk, 0 ],
					                       block_lut[ i_run, i_blk, 0 ] + block_len_vol
					                     ).astype( "int" )

					blk_data = np.mean( run_vtc[ blk_range ] )

					vox_data[ int( block_lut[ i_run, i_blk, 1 ] ) ].append( blk_data )

			stat[ i_voxel, : ] = scipy.stats.ttest_ind( vox_data[ 0 ],
			                                            vox_data[ 1 ]
			                                          )

		np.save( "%s-%s.npy" % ( paths[ "analysis" ][ "loc_stat" ],
		                         roi_name
		                       ),
		         arr = stat
		       )


def form_vtcs( paths, ana_conf, subj_conf, exp_conf, acq_conf ):
	"""Extracts the voxel time courses for each voxel in each ROI

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'glass_coherence.config'.
	ana_conf : dict
		Analysis configuration, as returned by 'get_analysis_conf' in
		'glass_coherence.config'.
	subj_conf : dict
		Subject configuration, as returned by 'get_subj_conf' in
		'glass_coherence.config', for this subject.
	exp_conf : dict
		Experiment configuration, as returned by 'get_exp_conf' in
		'glass_coherence.config'.

	"""

	for roi_name in ana_conf[ "rois" ]:

		# load the coordinates
		coords = np.load( "".join( [ paths[ "analysis" ][ "coords" ],
		                             "-", roi_name, ".npy"
		                           ]
		                         )
		                )

		# do a quick check to make sure the coords have the expected shape
		assert( coords.shape[ 0 ] == 3 )

		n_voxels = coords.shape[ 1 ]

		n_vols_per_run = int( exp_conf[ "n_blocks" ] *
		                      exp_conf[ "block_len_s" ] /
		                      acq_conf[ "tr_s" ]
		                    )

		vol_range_start = int( exp_conf[ "rej_start_blocks" ] *
		                       exp_conf[ "block_len_s" ] /
		                       acq_conf[ "tr_s" ]
		                     )

		vol_range_end = int( n_vols_per_run -
		                     ( exp_conf[ "rej_end_blocks" ] *
		                       exp_conf[ "block_len_s" ] /
		                       acq_conf[ "tr_s" ]
		                     )
		                   )

		vol_range = np.arange( vol_range_start, vol_range_end )

		vol_range += exp_conf[ "hrf_corr_vol" ]

		# initialise the vtc
		vtc = np.empty( ( n_vols_per_run,
		                  subj_conf[ "n_runs" ],
		                  n_voxels
		                )
		              )

		# fill with NaNs, to be safe
		vtc.fill( np.NAN )

		# loop through each unwarped image (run) file
		for ( i_run, run_path ) in enumerate( paths[ "func" ][ "uw" ] ):

			# load the run file
			run_img = nipy.load_image( "".join( [ run_path, ".nii" ] ) ).get_data()

			# iterate through each voxel in the roi
			for i_voxel in xrange( n_voxels ):

				# extract the voxel data (timecourse) at the voxel coordinate
				vox_data = run_img[ coords[ 0, i_voxel ],
				                    coords[ 1, i_voxel ],
				                    coords[ 2, i_voxel ],
				                    :
				                  ]

				# store the voxel data
				vtc[ :, i_run, i_voxel ] = vox_data


		# discard the unwanted volumes
		vtc = vtc[ vol_range, :, : ]

		# save the vtc
		np.save( "".join( [ paths[ "analysis" ][ "vtc" ],
		                   "-", roi_name
		                  ]
		                ),
		         arr = vtc
		       )

		# initialise the localiser vtc
		loc_vtc = np.empty( ( n_vols_per_run,
		                      subj_conf[ "n_loc_runs" ],
		                      n_voxels
		                    )
		                  )

		# fill with NaNs, to be safe
		loc_vtc.fill( np.NAN )

		# loop through each unwarped image (run) file
		for ( i_run, run_path ) in enumerate( paths[ "loc" ][ "uw" ] ):

			# load the run file
			run_img = nipy.load_image( "".join( [ run_path, ".nii" ] ) ).get_data()

			# iterate through each voxel in the roi
			for i_voxel in xrange( n_voxels ):

				# extract the voxel data (timecourse) at the voxel coordinate
				vox_data = run_img[ coords[ 0, i_voxel ],
				                    coords[ 1, i_voxel ],
				                    coords[ 2, i_voxel ],
				                    :
				                  ]

				# store the voxel data
				loc_vtc[ :, i_run, i_voxel ] = vox_data


		# discard the unwanted volumes
		loc_vtc = loc_vtc[ vol_range, :, : ]

		# save the vtc
		np.save( "".join( [ paths[ "analysis" ][ "loc_vtc" ],
		                   "-", roi_name
		                  ]
		                ),
		         arr = loc_vtc
		       )


def get_design( paths, acq_conf, subj_conf, exp_conf ):
	"""Extracts and writes the design matrix from the session log.

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'glass_coherence.config'.
	acq_conf : dict
		Acquisition configuration, as returned by 'get_acq_conf' in
		'glass_coherence.config'.
	subj_conf : dict
		Subject configuration, as returned by 'get_subj_conf' in
		'glass_coherence.config', for this subject.
	exp_conf : dict
		Experiment configuration, as returned by 'get_exp_conf' in
		'glass_coherence.config'.

	"""

	sess_log = scipy.io.loadmat( paths[ "design" ][ "log" ],
	                             squeeze_me = True,
	                             struct_as_record = False
	                           )

	# init the design matrix
	design = np.zeros( ( exp_conf[ "n_vols_per_run" ],
	                     exp_conf[ "n_cond" ],
	                     subj_conf[ "n_runs" ]
	                   )
	                 )

	for i_run in xrange( subj_conf[ "n_runs" ] ):

		evt = sess_log[ "runData" ][ i_run ].exp

		evt_lut = evt.evtLut

		# minus one to convert from Matlab's 1-based to Python's 0-based indices
		evt_i_cond = evt.iSeq - 1
		evt_i_onset_s = evt.iOnset - 1

		# cull the zero-contrast events
		valid_evt_lut = evt_lut[ evt_lut[ :, evt_i_cond ] > 1, : ]

		# extract the important info
		evt_cond = valid_evt_lut[ :, evt_i_cond ]
		evt_onset_s = valid_evt_lut[ :, evt_i_onset_s ]

		# convert the onset time to volumes
		evt_onset_vol = evt_onset_s / acq_conf[ "tr_s" ]

		# convert the condition to an index
		# -1 to get from [ 2, 5 ] to [ 1, 4 ], -1 to get to [ 0, 4 ]
		evt_cond = evt_cond - 1 - 1

		design[ evt_onset_vol.astype( 'int' ), evt_cond.astype( 'int' ), i_run ] = 1

	# discard the unwanted volume info
	design = design[ exp_conf[ "n_discard_vols" ]:, :, : ]

	assert( design.shape[ 0 ] == exp_conf[ "n_valid_vols_per_run" ] )

	for i_run in xrange( subj_conf[ "n_runs" ] ):

		shift_amount = subj_conf[ "onsets_adjust" ][ i_run ]

		design[ :, :, i_run ] = np.roll( design[ :, :, i_run ],
		                                 shift_amount,
		                                 axis = 0
		                               )

	# save the design matrix
	np.save( paths[ "design" ][ "matrix" ], arr = design )


def get_task( paths, exp_conf, subj_conf, ana_conf, acq_conf ):
	"""Extracts and writes the design matrix from the session log.

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'glass_coherence.config'.
	exp_conf : dict
		Experiment configuration, as returned by 'get_exp_conf' in
		'glass_coherence.config'.
	subj_conf : dict
		Subject configuration, as returned by 'get_subj_conf' in
		'glass_coherence.config', for this subject.


	!!UNFINISHED!!

	"""

	# calculate how many samples of task performance to take per run
	n_per_run = ( exp_conf[ "n_vols_per_run" ] *
	              acq_conf[ "tr_s" ] /
	              ana_conf[ "task_samp_rate_s" ]
	            )

	# initialise the task log
	# -first dimension is time, indexed at the task sampling rate
	# -second dimension is ( target present, response present, stimulus
	# condition )
	# -third dimension is acquistion run
	task_log = np.zeros( ( n_per_run,
	                       3,
	                       subj_conf[ "n_runs" ]
	                     )
	                   )

	# initialise the timer
	run_time = np.linspace( 0,
	                        exp_conf[ "n_vols_per_run" ] * acq_conf[ "tr_s" ],
	                        num = n_per_run,
	                        endpoint = False
	                      )

	# load the logfile for the acquisition session
	sess_log = scipy.io.loadmat( paths[ "design" ][ "log" ],
	                             squeeze_me = True,
	                             struct_as_record = False
	                           )

	for i_run in xrange( subj_conf[ "n_runs" ] ):

		# extract the log for this run
		run_log = sess_log[ "runData" ][ i_run ]

		# first, we want to build up our ingredients

		# -- 1. the target onset times
		# extract the task lut for this run
		task_lut = run_log.task.rsvp.taskLut
		# extract the time, in seconds, that the target events occurred
		target_s = task_lut[ task_lut[ :, -1 ] > 0,  # last dim indicates targets
		                     0  # first dim indicates time (s)
		                   ]

		# -- 2. the task response times
		resp_lut = run_log.task.resp
		# response vector was filled with NaNs to avoid 'growing' at runtime, so
		# only grab actual responses
		resp_s = resp_lut[ 1, ~np.isnan( resp_lut[ 0, : ] ) ]


		# -- 3. the active stimulus condition at each time
		stim_lut = run_log.exp.evtLut
		# stimulus time
		stim_s = stim_lut[ :, 1 ]
		# stimulus condition at each time
		stim_cond = stim_lut[ :, 0 ]

		# righto, now we have our ingredients we can build the vectors

		# iterate through the run at the desired task sample rate
		for ( i_time, evt_time ) in enumerate( run_time ):

			# time the event starts
			evt_s = evt_time
			# ... and ends
			evt_e = evt_s + ana_conf[ "task_samp_rate_s" ]

			# Q1 : was there a target on at this time?
			target_on = np.any( np.logical_and( target_s >= evt_s,
			                                    target_s < evt_e
			                                  )
			                  )

			# Q2 : was there a response at this time?
			resp_on = np.any( np.logical_and( resp_s >= evt_s,
			                                  resp_s < evt_e
			                                )
			                )

			# Q3 : what was the stimulus at this time?

			# work out which was the closest presentation in time (preceding)
			i_stim_evt = np.where( ( evt_s - stim_s ) >= 0 )[ 0 ][ -1 ]
			# get the condition for this event (minus one due to one-indexing)
			evt_cond = stim_cond[ i_stim_evt ] - 1

			# now we can update our container
			task_log[ i_time, :, i_run ] = [ target_on,
			                                 resp_on,
			                                 evt_cond
			                               ]

	np.save( paths[ "task" ][ "task_info" ],
	         task_log.astype( 'int' )
	       )
