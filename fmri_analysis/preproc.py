"""
Set of routines to pre-process the fMRI data for the natural scenes aperture
fMRI experiment.
"""

from __future__ import division

import os.path

import nipy
import numpy as np

import fmri_tools.preproc, fmri_tools.utils
import ns_aperture.fmri.loc


def convert( paths ):
	"""Converts the functionals and fieldmaps from dicom to nifti"""

	# aggregate the dicom directories
	raw_dirs = ( paths[ "func" ][ "raw_dirs" ] +
	             paths[ "loc" ][ "raw_dirs" ] +
	             paths[ "fmap" ][ "raw_mag_dirs" ] +
	             paths[ "fmap" ][ "raw_ph_dirs" ]
	           )

	# aggregate the output directories
	nii_dirs = ( paths[ "func" ][ "run_dirs" ] +
	             paths[ "loc" ][ "run_dirs" ] +
	             paths[ "fmap" ][ "fmap_dirs" ] +
	             paths[ "fmap" ][ "fmap_dirs" ]
	           )

	# aggregate the images paths
	img_paths = ( paths[ "func" ][ "orig_files" ] +
	              paths[ "loc" ][ "orig_files" ] +
	              paths[ "fmap" ][ "mag_files" ] +
	              paths[ "fmap" ][ "ph_files" ]
	            )

	# pull out the filenames
	img_names = [ os.path.split( img_path )[ 1 ]
	              for img_path in img_paths
	            ]

	# do the DCM -> NII conversion
	map( fmri_tools.preproc.dcm_to_nii,
	     raw_dirs,
	     nii_dirs,
	     img_names
	   )

	# generate the full paths (with assumed extension) of the newly-created nifti
	# files
	full_img_paths = [ "%s.nii" % img_path for img_path in img_paths ]

	# check that they are all unique
	assert( fmri_tools.utils.files_are_unique( full_img_paths ) )

	# files to go into the summary
	summ_paths = ( paths[ "func" ][ "orig_files" ] +
	               paths[ "loc" ][ "orig_files" ]
	             )

	# make a summary image from the files
	fmri_tools.preproc.gen_sess_summ_img( summ_paths,
	                                      paths[ "summ" ][ "orig_summ_file" ]
	                                    )


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
	orig_paths = [ paths[ im_type ][ "orig_files" ][ i_run - 1 ]
	               for i_run, im_type in run_order
	             ]
	corr_paths = [ paths[ im_type ][ "corr_files" ][ i_run - 1 ]
	               for i_run, im_type in run_order
	             ]

	# pull out the important information from the config
	slice_order = conf[ "acq" ][ "slice_order" ]
	tr_s = conf[ "acq" ][ "tr_s" ]
	slice_info = ( conf[ "acq" ][ "slice_axis" ],
	               conf[ "acq" ][ "slice_acq_dir" ]
	             )
	st_correct = conf[ "preproc" ][ "st_correct" ]

	# run the motion correction algorithm (slow)
	motion_est = fmri_tools.preproc.correct_st_motion( orig_paths,
	                                                   corr_paths,
	                                                   slice_order,
	                                                   tr_s,
	                                                   slice_info,
	                                                   st_correct = st_correct
	                                                 )

	# save the estimated motion parameters
	np.save( paths[ "summ" ][ "mot_est_file" ],
	         arr = motion_est
	       )

	# make a summary image from the corrected files
	fmri_tools.preproc.gen_sess_summ_img( corr_paths,
	                                      paths[ "summ" ][ "corr_summ_file" ]
	                                    )


def fieldmaps( paths, conf, subj_conf ):
	"""Prepare the fieldmaps"""

	# duplicate the delta TE for each fieldmap acquired
	delta_te_ms = ( [ conf[ "acq" ][ "delta_te_ms" ] ] *
	                subj_conf[ "n_fmaps" ]
	              )

	map( fmri_tools.preproc.make_fieldmap,
	     paths[ "fmap" ][ "mag_files" ],
	     paths[ "fmap" ][ "ph_files" ],
	     paths[ "fmap" ][ "fmap_files" ],
	     delta_te_ms
	   )


def unwarp( paths, conf ):
	"""Uses the fieldmaps to unwarp the functional images and create a mean image
	of all the unwarped functional images.
	"""

	# combine the experiment and localiser functional info
	func_fmap = paths[ "func" ][ "fmap_files" ] + paths[ "loc" ][ "fmap_files" ]
	func_corr = paths[ "func" ][ "corr_files" ] + paths[ "loc" ][ "corr_files" ]
	func_uw = paths[ "func" ][ "uw_files" ] + paths[ "loc" ][ "uw_files" ]

	# duplicate the dwell time and phase encode direction for each image
	dwell_ms = [ conf[ "acq" ][ "dwell_ms" ] ] * len( func_corr )
	ph_encode_dir = [ conf[ "acq" ][ "ph_encode_dir" ] ] * len( func_corr )

	interp = [ "spline" ] * len( func_corr )

	# perform the unwarping
	map( fmri_tools.preproc.unwarp,
	     func_corr,
	     func_fmap,
	     func_uw,
	     dwell_ms,
	     ph_encode_dir,
	     interp
	   )

	# create a mean image of the unwarped data
	fmri_tools.preproc.mean_image( func_uw,
	                               paths[ "summ" ][ "mean_file" ]
	                             )

	# produce a summary image
	fmri_tools.preproc.gen_sess_summ_img( func_corr,
	                                      paths[ "summ" ][ "uw_summ_file" ]
	                                    )


def make_roi_images( paths, conf ):
	"""Converts the ROI matlab files to nifti images in register with the
	subject's anatomical.
	"""

	n_rois = len( conf[ "ana" ][ "rois" ] )

	# python nifti loader needs the extension
	anat_path = [ "%s.nii" % paths[ "anat" ][ "anat_file" ] ] * n_rois

	roi_paths = [ "%s.nii" % roi_path
	              for roi_path in paths[ "roi" ][ "orig_files" ]
	            ]

	roi_ax = [ conf[ "ana" ][ "roi_ax" ] ] * n_rois
	roi_ax_order = [ conf[ "ana" ][ "roi_ax_order" ] ] * n_rois

	# convert the ROI mat files to nifti images
	map( fmri_tools.preproc.roi_to_nii,
	     paths[ "roi" ][ "mat_files" ],
	     anat_path,
	     roi_paths,
	     roi_ax,
	     roi_ax_order
	   )


def prepare_rois( paths, conf ):
	"""Extracts and writes the coordinates of the gray matter and each ROI"""

	# load the gray matter mask
	gray_img = nipy.load_image( "%s.nii" % paths[ "anat" ][ "gray_rs_file" ]
	                          ).get_data()

	# get the coordinates, as 3 x N
	gray_coords = np.array( np.nonzero( gray_img ) )

	# and save
	np.save( paths[ "roi" ][ "gray_coord_file" ],
	         arr = gray_coords
	       )

	n_gray_voxels = gray_coords.shape[ 1 ]

	# this image will hold the index, in gray_coords, of each gray matter voxel
	# over the volume - this will be used as a lookup table for the ROIs
	i_gray_img = np.empty( gray_img.shape )
	i_gray_img.fill( np.NAN )

	for i_coord in xrange( n_gray_voxels ):

		i_gray_img[ gray_coords[ 0, i_coord ],
		            gray_coords[ 1, i_coord ],
		            gray_coords[ 2, i_coord ]
		          ] = i_coord


	for i_roi in xrange( len( conf[ "ana" ][ "rois" ] ) ):

		roi = nipy.load_image( "%s.nii" % paths[ "roi" ][ "rs_files" ][ i_roi ]
		                     ).get_data()

		roi_coords = np.array( np.nonzero( roi ) )

		n_roi_coords = roi_coords.shape[ 1 ]

		roi_gray_coords = np.empty( ( n_roi_coords ) )
		roi_gray_coords.fill( np.NAN )

		for i_coord in xrange( n_roi_coords ):

			# find the gray coordinate index for this ROI voxel
			roi_gray_coords[ i_coord ] = i_gray_img[ roi_coords[ 0, i_coord ],
			                                         roi_coords[ 1, i_coord ],
			                                         roi_coords[ 2, i_coord ]
			                                       ]

		i_not_in_gray = np.where( np.isnan( roi_gray_coords ) )[ 0 ]

		roi_gray_coords = np.delete( roi_gray_coords, i_not_in_gray )

		missing = len( i_not_in_gray ) / n_roi_coords

		if missing > 0.01:
			print "Warning: more than 1% of ROI coords not within gray matter mask"

		# and save
		np.save( paths[ "roi" ][ "coord_files" ][ i_roi ],
		         arr = roi_gray_coords.astype( "int" )
		       )


def form_vtcs( paths, conf, subj_conf ):
	"""Extracts the voxel time courses for each voxel in the gray matter"""

	gray_coords = np.load( paths[ "roi" ][ "gray_coord_file" ] )

	vtc = np.empty( ( conf[ "exp" ][ "n_vols_per_run" ],
	                  subj_conf[ "n_runs" ],
	                  gray_coords.shape[ 1 ]
	                )
	              )
	vtc.fill( np.NAN )

	# loop through each unwarped image (run) file
	for ( i_run, run_path ) in enumerate( paths[ "func" ][ "uw_files" ] ):

		# load the run data from the image file file
		run_img = nipy.load_image( "%s.nii" % run_path ).get_data()

		for i_coord in xrange( gray_coords.shape[ 1 ] ):

			run_vtc = run_img[ gray_coords[ 0, i_coord ],
			                   gray_coords[ 1, i_coord ],
			                   gray_coords[ 2, i_coord ],
			                   :
			                 ]

			run_vtc = fmri_tools.preproc.hp_filter( run_vtc,
			                                        filt_type = "fir",
			                                        as_psc = False,
			                                        rem_mean = False,
			                                        tr_s = conf[ "acq" ][ "tr_s" ],
			                                      )[ 0 ]

			vtc[ :, i_run, i_coord ] = run_vtc

	# save the vtc
	np.save( "%s-gray.npy" % paths[ "ana_exp" ][ "vtc_file" ],
	         arr = vtc
	       )

	# LOCALISERS
	loc_vtc = np.empty( ( conf[ "exp" ][ "loc_n_vols_per_run" ],
	                      subj_conf[ "n_loc_runs" ],
	                      gray_coords.shape[ 1 ]
	                    )
	                  )
	loc_vtc.fill( np.NAN )

	# loop through each unwarped image (run) file
	for ( i_run, run_path ) in enumerate( paths[ "loc" ][ "uw_files" ] ):

		# load the run data from the image file file
		run_img = nipy.load_image( "%s.nii" % run_path ).get_data()

		for i_coord in xrange( gray_coords.shape[ 1 ] ):

			run_vtc = run_img[ gray_coords[ 0, i_coord ],
			                   gray_coords[ 1, i_coord ],
			                   gray_coords[ 2, i_coord ],
			                   :
			                 ]

			run_vtc = fmri_tools.preproc.hp_filter( run_vtc,
			                                        filt_type = "fir",
			                                        as_psc = False,
			                                        rem_mean = False,
			                                        tr_s = conf[ "acq" ][ "tr_s" ],
			                                      )[ 0 ]

			loc_vtc[ :, i_run, i_coord ] = run_vtc

	# save the vtc
	np.save( "%s-gray.npy" % paths[ "ana_loc" ][ "vtc_file" ],
	         arr = loc_vtc
	       )


def cull_voxels( paths, conf, subj_conf ):
	"""Cull any bad voxels"""

	vtc = np.load( "%s-gray.npy" % paths[ "ana_exp" ][ "vtc_file" ] )
	loc_vtc = np.load( "%s-gray.npy" % paths[ "ana_loc" ][ "vtc_file" ] )

	i_good = np.ones( ( vtc.shape[ -1 ] ) ).astype( "bool" )

	for i_voxel in xrange( len( i_good ) ):

		if np.any( vtc[ :, :, i_voxel ] <= 0 ):
			i_good[ i_voxel ] = False

		if np.any( loc_vtc[ :, :, i_voxel ] <= 0 ):
			i_good[ i_voxel ] = False

	vtc = vtc[ :, :, i_good ]
	np.save( "%s-gray.npy" % paths[ "ana_exp" ][ "vtc_file" ], vtc )

	loc_vtc = loc_vtc[ :, :, i_good ]
	np.save( "%s-gray.npy" % paths[ "ana_loc" ][ "vtc_file" ], loc_vtc )

	# cull from the coordinate list
	gray_coords = np.load( paths[ "roi" ][ "gray_coord_file" ] )

	gray_coords = gray_coords[ :, i_good ]

	np.save( paths[ "roi" ][ "gray_coord_file" ], gray_coords )

	# ... and from the ROIs
	for i_roi in xrange( len( conf[ "ana" ][ "rois" ] ) ):

		roi_coords = np.load( paths[ "roi" ][ "coord_files" ][ i_roi ] )

		i_roi_good = i_good[ roi_coords ]

		roi_coords = roi_coords[ i_roi_good ]

		np.save( paths[ "roi" ][ "coord_files" ][ i_roi ], roi_coords )


def get_design( paths, conf, subj_conf ):
	"""Extracts and writes the design matrix from the session log"""

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

		log_path = "%s%d.npy" % ( paths[ "log" ][ "seq_base" ],
		                           i_run + 1
		                         )

		log = np.load( log_path )

		# iterate through only the block indices we care about
		for ( i_block, block) in enumerate( block_range ):

			# find the first event corresponding to this block indice
			# because they're indices, need to add 1 to match the 1-based storage in
			# the log
			i_block_start_evt = np.where( log[ :, 1 ] == ( block + 1 ) )[ 0 ][ 0 ]

			# find the start volume for this block
			block_start_vol = log[ i_block_start_evt, 0 ] / conf[ "acq" ][ "tr_s" ]

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

	np.save( paths[ "log" ][ "design" ], design )


	# localiser
	loc_ord = ( "AB", "BA" )

	loc_design = np.empty( ( conf[ "exp" ][ "loc_n_valid_blocks" ],
	                         subj_conf[ "n_loc_runs" ],
	                         2
	                       )
	                     )
	loc_design.fill( np.NAN )

	loc_block_range = np.arange( 0, conf[ "exp" ][ "loc_n_valid_blocks" ] )

	for ( i_run, run_loc_ord ) in enumerate( loc_ord ):

		loc_seq = ns_aperture.fmri.loc.get_seq( conf,
		                                        run_loc_ord
		                                      )

		for ( i_block, block ) in enumerate( loc_block_range ):

			# find the first event corresponding to this block indice
			# because they're indices, need to add 1 to match the 1-based storage in
			# the log
			i_block_start_evt = np.where( loc_seq[ :, 1 ] == ( block + 1 ) )[ 0 ][ 0 ]

			# find the start volume for this block
			block_start_vol = loc_seq[ i_block_start_evt, 0 ] / conf[ "acq" ][ "tr_s" ]

			# pull out the condition type for the block
			block_cond = loc_seq[ i_block_start_evt, 2 ]

			# store
			loc_design[ i_block, i_run, 0 ] = block_start_vol
			loc_design[ i_block, i_run, 1 ] = block_cond

	# make sure there are equal numbers of each condition type, as expected
	assert( np.all( np.sum( loc_design[ :, :, 1 ] == 0 ) ==
	                np.sum( loc_design[ :, :, 1 ] == 1 )
	              )
	      )

	np.save( paths[ "log" ][ "loc_design" ], loc_design )
