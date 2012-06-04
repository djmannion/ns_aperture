"""
Set of routines to pre-process the fMRI data for the natural scenes aperture
fMRI experiment.
"""

from __future__ import division

import os.path

import nipy
import numpy as np
import scipy.io

import fmri_tools.preproc, fmri_tools.utils

def convert( paths ):
	"""Converts the functionals and fieldmaps from dicom to nifti"""

	# aggregate the dicom directories
	dcm_dirs = ( paths[ "func" ][ "dcm_dirs" ] +
	             paths[ "loc" ][ "dcm_dirs" ] +
	             paths[ "fmap" ][ "dcm_mag_dirs" ] +
	             paths[ "fmap" ][ "dcm_ph_dirs" ]
	           )

	# aggregate the output directories
	nii_dirs = ( paths[ "func" ][ "dirs" ] +
	             paths[ "loc" ][ "dirs" ] +
	             paths[ "fmap" ][ "dirs" ] +
	             paths[ "fmap" ][ "dirs" ]
	           )

	# aggregate the images paths
	img_paths = ( paths[ "func" ][ "orig" ] +
	              paths[ "loc" ][ "orig" ] +
	              paths[ "fmap" ][ "mag" ] +
	              paths[ "fmap" ][ "ph" ]
	            )

	# pull out the filenames
	img_names = [ os.path.split( img_path )[ 1 ]
	              for img_path in img_paths
	            ]

	# do the DCM -> NII conversion
	map( fmri_tools.preproc.dcm_to_nii,
	     dcm_dirs,
	     nii_dirs,
	     img_names
	   )

	# generate the full paths (with assumed extension) of the newly-created nifti
	# files
	full_img_paths = [ "%s.nii" % img_path for img_path in img_paths ]

	# check that they are all unique
	assert( fmri_tools.utils.files_are_unique( full_img_paths ) )


def st_motion_correct( paths, conf, subj_conf ):
	"""Performs slice-timing and motion correction"""

	# get the order of runs to pass to the correction algorithm
	# this is done because the algorithm realigns all to the first entry, which
	# normally corresponds to the first run acquired, but we might want them to
	# be aligned with a different run - one closer to a fieldmap, for example
	run_order = subj_conf[ "run_st_mot_order" ]

	# reorder the paths
	# (the -1 is because the runs are specified in subj_conf in a one-based
	# index; ie. run 1 is the first run)
	orig_paths = [ paths[ im_type ][ "orig" ][ i_run - 1 ]
	               for i_run, im_type in run_order
	             ]
	corr_paths = [ paths[ im_type ][ "corr" ][ i_run - 1 ]
	               for i_run, im_type in run_order
	             ]

	# pull out the important information from the config
	slice_order = conf[ "acq" ][ "slice_order" ]
	tr_s = conf[ "acq" ][ "tr_s" ]
	slice_info = ( conf[ "acq" ][ "slice_axis" ],
	               conf[ "acq" ][ "slice_acq_dir" ]
	             )

	# run the motion correction algorithm (slow)
	motion_est = fmri_tools.preproc.correct_st_motion( orig_paths,
	                                                   corr_paths,
	                                                   slice_order,
	                                                   tr_s,
	                                                   slice_info,
	                                                 )

	# save the estimated motion parameters
	np.save( paths[ "func" ][ "motion_estimates" ],
	         arr = motion_est
	       )

	# make a summary image from the corrected files
	fmri_tools.preproc.gen_sess_summ_img( corr_paths,
	                                      paths[ "func" ][ "corr_summ" ]
	                                    )


def fieldmaps( paths, conf, subj_conf ):
	"""Prepare the fieldmaps"""

	# duplicate the delta TE for each fieldmap acquired
	delta_te_ms = ( [ conf[ "acq" ][ "delta_te_ms" ] ] *
	                subj_conf[ "n_fmaps" ]
	              )

	map( fmri_tools.preproc.make_fieldmap,
	     paths[ "fmap" ][ "mag" ],
	     paths[ "fmap" ][ "ph" ],
	     paths[ "fmap" ][ "fmap" ],
	     delta_te_ms
	   )


def unwarp( paths, conf ):
	"""Uses the fieldmaps to unwarp the functional images and create a mean image
	of all the unwarped functional images.
	"""

	# combine the experiment and localiser functional info
	func_corr = paths[ "func" ][ "corr" ] + paths[ "loc" ][ "corr" ]
	func_fmap = paths[ "func" ][ "fmap" ] + paths[ "loc" ][ "fmap" ]
	func_uw = paths[ "func" ][ "uw" ] + paths[ "loc" ][ "uw" ]

	# duplicate the dwell time and phase encode direction for each image
	dwell_ms = [ conf[ "acq" ][ "dwell_ms" ] ] * len( func_corr )
	ph_encode_dir = [ conf[ "acq" ][ "ph_encode_dir" ] ] * len( func_corr )

	# perform the unwarping
	map( fmri_tools.preproc.unwarp,
	     func_corr,
	     func_fmap,
	     func_uw,
	     dwell_ms,
	     ph_encode_dir
	   )

	# create a mean image of the unwarped data
	fmri_tools.preproc.mean_image( func_uw,
	                               paths[ "func" ][ "mean" ]
	                             )

	# produce a summary image
	fmri_tools.preproc.gen_sess_summ_img( func_uw,
	                                      paths[ "func" ][ "uw_summ" ]
	                                    )


def make_roi_images( paths, conf ):
	"""Converts the ROI matlab files to nifti images in register with the
	subject's anatomical.
	"""

	n_rois = len( conf[ "ana" ][ "rois" ] )

	# python nifti loader needs the extension
	anat_path = [ "%s.nii" % paths[ "anat" ][ "anat" ] ] * n_rois

	roi_paths = [ "%s.nii" % roi_path
	              for roi_path in paths[ "roi" ][ "orig" ]
	            ]

	roi_ax = [ conf[ "ana" ][ "roi_ax" ] ] * n_rois
	roi_ax_order = [ conf[ "ana" ][ "roi_ax_order" ] ] * n_rois

	# convert the ROI mat files to nifti images
	map( fmri_tools.preproc.roi_to_nii,
	     paths[ "roi" ][ "mat" ],
	     anat_path,
	     roi_paths,
	     roi_ax,
	     roi_ax_order
	   )


def prepare_rois( paths, conf ):
	"""Extracts and writes the coordinates of each ROI"""

	for ( i_roi, roi_name ) in enumerate( conf[ "ana" ][ "rois" ] ):

		# load the ROI image
		roi = nipy.load_image( "%s.nii" % paths[ "roi" ][ "rs" ][ i_roi ]
		                     ).get_data()

		# get rid of any NaNs
		roi[ np.isnan( roi) ] = 0

		# extract the coordinate list
		coords = np.nonzero( roi )

		# and save
		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "coords" ], roi_name ),
		         arr = coords
		       )


def form_vtcs( paths, conf, subj_conf ):
	"""Extracts the voxel time courses for each voxel in each ROI"""

	for roi_name in conf[ "ana" ][ "rois" ]:

		# load the coordinates for the roi
		coords = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "coords" ], roi_name ) )

		n_voxels = coords.shape[ 1 ]

		# initialise the vtc
		vtc = np.empty( ( conf[ "exp" ][ "n_vols_per_run" ],
		                  subj_conf[ "n_runs" ],
		                  n_voxels
		                )
		              )

		# fill with NaNs, to be safe
		vtc.fill( np.NAN )

		# loop through each unwarped image (run) file
		for ( i_run, run_path ) in enumerate( paths[ "func" ][ "uw" ] ):

			# load the run file
			run_img = nipy.load_image( "%s.nii" % run_path ).get_data()

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
		vtc = vtc[ conf[ "exp" ][ "run_range" ], :, : ]

		# check that it has been filled up correctly
		assert( np.sum( np.isnan( vtc ) ) == 0 )

		# save the vtc
		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "vtc" ], roi_name ),
		         arr = vtc
		       )

		# initialise the localiser vtc
		loc_vtc = np.empty( ( conf[ "exp" ][ "n_vols_per_run" ],
		                      subj_conf[ "n_loc_runs" ],
		                      n_voxels
		                    )
		                  )

		# fill with NaNs, to be safe
		loc_vtc.fill( np.NAN )

		# loop through each unwarped image (run) file
		for ( i_run, run_path ) in enumerate( paths[ "loc" ][ "uw" ] ):

			# load the run file
			run_img = nipy.load_image( "%s.nii" % run_path ).get_data()

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

		# discard the unwanted volumes and compensate for the HRF delay
		loc_vtc = loc_vtc[ conf[ "exp" ][ "run_range" ], :, : ]

		# check that it has been filled up correctly
		assert( np.sum( np.isnan( loc_vtc ) ) == 0 )

		# save the localiser vtc
		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "loc_vtc" ], roi_name ),
		         arr = loc_vtc
		       )


def get_design( paths, conf, subj_conf ):
	"""Extracts and writes the design matrix from the session log.

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'ns_aperture.config'.
	conf : dict
		Experiment configuration, as returned by 'get_conf' in
		'ns_aperture.config'.
	subj_conf : dict
		Subject configuration, as returned by 'get_subj_conf' in
		'ns_aperture.config', for this subject.

	"""

	design = np.empty( ( conf[ "exp" ][ "n_valid_blocks" ],
	                     subj_conf[ "n_runs" ],
	                     2
	                   )
	                 )
	design.fill( np.NAN )

	# carries the block indices (0-based) we care about
	block_range = np.arange( conf[ "exp" ][ "rej_start_blocks" ],
	                         conf[ "exp" ][ "n_blocks" ] -
	                         conf[ "exp" ][ "rej_end_blocks" ]
	                       )

	for i_run in xrange( subj_conf[ "n_runs" ] ):

		log_path = os.path.join( paths[ "design" ][ "log_dir" ],
		                         "%s_ns_aperture_fmri_seq_%d.npy" % (
		                         subj_conf[ "subj_id" ],
		                         i_run + 1
		                         )
		                       )

		log = np.load( log_path )

		# iterate through only the block indices we care about
		for ( i_block, block) in enumerate( block_range ):

			# find the first event corresponding to this block indice
			# because they're indices, need to add 1 to match the 1-based storage in
			# the log
			i_block_start_evt = np.where( log[ :, 1 ] == ( block + 1 ) )[ 0 ][ 0 ]

			# find the start volume for this block, discounting invalid blocks
			block_start_vol = i_block * conf[ "exp" ][ "n_vols_per_blk" ]

			# pull out the condition type for the block
			block_cond = log[ i_block_start_evt, 2 ]

			# store
			design[ i_block, i_run, 0 ] = block_start_vol
			design[ i_block, i_run, 1 ] = block_cond

	# make sure there are equal numbers of each condition type, as expected
	assert( np.all( np.sum( design[ :, :, 1 ] == 0 ) ==
	                np.sum( design[ :, :, 1 ] == 1 )
	              )
	      )

	np.save( paths[ "design" ][ "design" ], design )

	# localiser follows the same main design, only for less runs
	loc_design = design[ :, :subj_conf[ "n_loc_runs" ], : ]

	np.save( paths[ "design" ][ "loc_design" ], loc_design )

	return design


def localiser_analysis( paths, conf ):
	"""Performs statistical analysis of the localiser data.

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'ns_aperture.config'.
	conf : dict
		Experiment configuration, as returned by 'get_conf' in
		'ns_aperture.config'.

	"""

	# load the design info for the localiser
	design = np.load( paths[ "design" ][ "loc_design" ] )

	for roi_name in conf[ "ana" ][ "rois" ]:

		# load the (trimmed and HRF corrected) vtc
		loc_vtc = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "loc_vtc" ],
		                                   roi_name
		                                 )
		                 )

		( _, n_runs, n_voxels ) = loc_vtc.shape

		stat = np.empty( ( n_voxels, 2 ) )
		stat.fill( np.NAN )

		for i_voxel in xrange( n_voxels ):

			vox_data = [ [], [] ]

			for i_run in xrange( n_runs ):

				run_vtc = loc_vtc[ :, i_run, i_voxel ]
				run_vtc = fmri_tools.preproc.hp_filter( run_vtc,
				                                        poly_ord = conf[ "ana" ][ "poly_ord" ]
				                                      )[ 0 ]

				run_design = design[ :, i_run, : ].astype( "int" )

				for i_blk in xrange( run_design.shape[ 0 ] ):

					blk_vol_range = np.arange( run_design[ i_blk, 0 ],
					                           run_design[ i_blk, 0 ] +
					                           conf[ "exp" ][ "n_vols_per_blk" ]
					                         ).astype( "int" )

					blk_data = np.mean( run_vtc[ blk_vol_range ] )

					vox_data[ run_design[ i_blk, 1 ] ].append( blk_data )

			stat[ i_voxel, : ] = scipy.stats.ttest_ind( vox_data[ 0 ],
			                                            vox_data[ 1 ]
			                                          )

		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "loc_stat" ],
		                         roi_name
		                       ),
		         arr = stat
		       )


def voxel_selection( paths, conf ):
	"""Finds the voxels in each ROI that are significantly activated by the
	   localiser, and extracts their coordinates and contrains the vtcs.

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'ns_aperture.config'.
	conf : dict
		Experiment configuration, as returned by 'get_conf' in
		'ns_aperture.config'.

	"""

	for roi_name in conf[ "ana" ][ "rois" ]:

		# all the voxels for this ROI
		full_coords = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "coords" ],
		                                       roi_name
		                                     )
		                     )

		# the localiser analysis
		loc_stat = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "loc_stat" ],
		                                    roi_name
		                                  )
		                  )

		# make sure they have the right number of voxels
		assert( full_coords.shape[ 1 ] == loc_stat.shape[ 0 ] )

		# activation probability for each voxel
		p_val = loc_stat[ :, 1 ]

		# perform the thresholding
		i_valid = ( p_val < conf[ "ana" ][ "loc_p_thresh" ] )

		# cull the coords not above threshold
		sel_coords = full_coords[ :, i_valid ]

		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "coords_sel" ], roi_name ),
		         sel_coords
		       )

		# apply culling to the vtcs
		vtc = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "vtc" ], roi_name ) )
		vtc = vtc[ :, :, i_valid ]
		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "vtc_sel" ], roi_name ),
		         vtc
		       )

		loc_vtc = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "loc_vtc" ], roi_name ) )
		loc_vtc = loc_vtc[ :, :, i_valid ]
		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "loc_vtc_sel" ], roi_name ),
		         loc_vtc
		       )


def avg_vtcs( paths, conf ):
	"""Averages the timecourses over all the (selected) voxels in a ROI and
	   performs high-pass filtering.

	Parameters
	----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'ns_aperture.config'.
	conf : dict
		Experiment configuration, as returned by 'get_analysis_conf' in
		'ns_aperture.config'.

	"""

	p_ord = conf[ "ana" ][ "poly_ord" ]

	for roi_name in conf[ "ana" ][ "rois" ]:

		# load the (post voxel selection) vtc
		vtc = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "vtc_sel" ],
		                               roi_name
		                             )
		             )

		# average over voxels
		vtc = np.mean( vtc, axis = 2 )

		filt_vtc = np.empty( vtc.shape )
		filt_vtc.fill( np.NAN )

		for i_run in xrange( filt_vtc.shape[ 1 ] ):

			filt_vtc[ :, i_run ] = fmri_tools.preproc.hp_filter( vtc[ :, i_run ],
			                                                     poly_ord = p_ord
			                                                   )[ 0 ]

		assert( not np.any( np.isnan( filt_vtc ) ) )

		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "vtc_avg" ], roi_name ),
		         filt_vtc
		       )
