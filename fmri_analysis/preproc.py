"""
Set of routines to pre-process the fMRI data for the natural scenes aperture
fMRI experiment.
"""

from __future__ import division

import os.path

import numpy as np

import fmri_tools.preproc, fmri_tools.utils


def convert( paths, conf ):
	"""Converts the functionals and fieldmaps from dicom to nifti"""

	# aggregate the dicom directories
	raw_dirs = ( paths[ "func_exp" ][ "raw_dirs" ] +
	             paths[ "func_loc" ][ "raw_dirs" ] +
	             paths[ "fmap" ][ "raw_mag_dirs" ] +
	             paths[ "fmap" ][ "raw_ph_dirs" ]
	           )

	# aggregate the output directories
	nii_dirs = ( paths[ "func_exp" ][ "run_dirs" ] +
	             paths[ "func_loc" ][ "run_dirs" ] +
	             paths[ "fmap" ][ "fmap_dirs" ] +
	             paths[ "fmap" ][ "fmap_dirs" ]
	           )

	# aggregate the images paths
	img_paths = ( paths[ "func_exp" ][ "orig_files" ] +
	              paths[ "func_loc" ][ "orig_files" ] +
	              paths[ "fmap" ][ "mag_files" ] +
	              paths[ "fmap" ][ "ph_files" ]
	            )

	# pull out the filenames
	img_names = [ os.path.split( img_path )[ 1 ]
	              for img_path in img_paths
	            ]

	for i_dir in xrange( len( raw_dirs ) ):

		fmri_tools.preproc.dcm_to_nii( raw_dirs[ i_dir ],
		                               nii_dirs[ i_dir ],
		                               img_names[ i_dir ],
		                               reorient_dim = conf[ "acq" ][ "ras" ],
		                               log_path = paths[ "summ" ][ "log_file" ]
		                             )


	# generate the full paths (with assumed extension) of the newly-created nifti
	# files
	full_img_paths = [ "%s.nii" % img_path for img_path in img_paths ]

	# check that they are all unique
	assert( fmri_tools.utils.files_are_unique( full_img_paths ) )

	# files to go into the summary
	summ_paths = ( paths[ "func_exp" ][ "orig_files" ] +
	               paths[ "func_loc" ][ "orig_files" ]
	             )

	# make a summary image from the files
	fmri_tools.preproc.gen_sess_summ_img( summ_paths,
	                                      paths[ "summ" ][ "orig_summ_file" ],
	                                      paths[ "summ" ][ "log_file" ]
	                                    )


def st_motion_correct( paths, conf ):
	"""Performs slice-timing and motion correction"""

	# get the order of runs to pass to the correction algorithm
	# this is done because the algorithm realigns all to the first entry, which
	# normally corresponds to the first run acquired, but we might want them to
	# be aligned with a different run - one closer to a fieldmap, for example
	run_order = conf[ "subj" ][ "run_st_mot_order" ]

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

	# run the motion correction algorithm (slow)
	motion_est = fmri_tools.preproc.correct_st_motion( orig_paths,
	                                                   corr_paths,
	                                                   slice_order,
	                                                   tr_s,
	                                                   slice_info
	                                                 )

	# save the estimated motion parameters
	np.save( paths[ "summ" ][ "mot_est_file" ],
	         arr = motion_est
	       )

	# make a summary image from the corrected files
	fmri_tools.preproc.gen_sess_summ_img( corr_paths,
	                                      paths[ "summ" ][ "corr_summ_file" ]
	                                    )


def fieldmaps( paths, conf ):
	"""Prepare the fieldmaps"""

	for i_fmap in xrange( conf[ "subj" ][ "n_fmaps" ] ):

		fmri_tools.preproc.make_fieldmap( paths[ "fmap" ][ "mag_files" ][ i_fmap ],
		                                  paths[ "fmap" ][ "ph_files" ][ i_fmap ],
		                                  paths[ "fmap" ][ "fmap_files" ][ i_fmap ],
		                                  conf[ "acq" ][ "delta_te_ms" ],
		                                  log_path = paths[ "summ" ][ "log_file" ]
		                                )


def unwarp( paths, conf ):
	"""Uses the fieldmaps to unwarp the functional images and create a mean image
	of all the unwarped functional images.
	"""

	# combine the experiment and localiser functional info
	func_fmap = ( paths[ "func_exp" ][ "fmap_files" ] +
	              paths[ "func_loc" ][ "fmap_files" ]
	            )

	func_corr = ( paths[ "func_exp" ][ "corr_files" ] +
	              paths[ "func_loc" ][ "corr_files" ]
	            )

	func_uw = ( paths[ "func_exp" ][ "uw_files" ] +
	            paths[ "func_loc" ][ "uw_files" ]
	          )

	for _ in xrange( conf[ "subj" ][ "n_runs" ] ):

		fmri_tools.preproc.unwarp( func_corr,
		                           func_fmap,
		                           func_uw,
		                           conf[ "acq" ][ "dwell_ms" ],
		                           conf[ "acq" ][ "ph_encode_dir" ],
		                           log_path = paths[ "summ" ][ "log_file" ]
		                         )

	# create a mean image of the unwarped data
	fmri_tools.preproc.mean_image( func_uw,
	                               paths[ "summ" ][ "mean_file" ],
	                               log_path = paths[ "summ" ][ "log_file" ]
	                             )

	# produce a summary image
	fmri_tools.preproc.gen_sess_summ_img( func_uw,
	                                      paths[ "summ" ][ "uw_summ_file" ],
	                                      log_path = paths[ "summ" ][ "log_file" ]
	                                    )
