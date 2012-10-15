"""
Set of routines to pre-process the fMRI data for the natural scenes aperture
fMRI experiment.
"""

from __future__ import division

import os.path

import fmri_tools.preproc, fmri_tools.utils


def convert( paths, conf ):
	"""Converts the functionals and fieldmaps from dicom to nifti"""

	raw_dirs = ( paths[ "func" ][ "raw_dirs" ] +
	             paths[ "fmap" ][ "raw_mag_dirs" ] +
	             paths[ "fmap" ][ "raw_ph_dirs" ]
	           )

	nii_files = ( paths[ "func" ][ "orig_files" ] +
	              paths[ "fmap" ][ "mag_files" ] +
	              paths[ "fmap" ][ "ph_files" ]
	            )

	for ( raw_dir, nii_file ) in zip( raw_dirs, nii_files ):

		fmri_tools.preproc.dcm_to_nii( raw_dir,
		                               "%s.nii" % nii_file,
		                               reorient_dim = conf[ "acq" ][ "ras" ],
		                               log_path = paths[ "summ" ][ "log_file" ]
		                             )

	# generate the full paths (with assumed extension) of the newly-created nifti
	# files
	full_nii_paths = [ "%s.nii" % nii_file for nii_file in nii_files ]

	# check that they are all unique
	assert( fmri_tools.utils.files_are_unique( full_nii_paths ) )

	# make a summary image from the files
	fmri_tools.preproc.gen_sess_summ_img( paths[ "func" ][ "orig_files" ],
	                                      paths[ "summ" ][ "orig_summ_file" ],
	                                      log_path = paths[ "summ" ][ "log_file" ]
	                                    )


def motion_correct( paths, conf ):
	"""Performs motion correction"""

	# the index to the run that will be the 'base', to be corrected to
	# it is stored in one-based, hence the minus 1 to get to an index
	i_mc_base = conf[ "subj" ][ "mot_base" ] - 1

	mc_base = "%s.nii[0]" % paths[ "func" ][ "orig_files" ][ i_mc_base ]

	fmri_tools.preproc.mot_correct( paths[ "func" ][ "orig_files" ],
	                                paths[ "func" ][ "corr_files" ],
	                                mc_base,
	                                mc_path = paths[ "summ" ][ "mot_est_file" ],
	                                log_path = paths[ "summ" ][ "log_file" ]
	                              )

	# make a summary image from the corrected files
	fmri_tools.preproc.gen_sess_summ_img( paths[ "func" ][ "corr_files" ],
	                                      paths[ "summ" ][ "corr_summ_file" ],
	                                      log_path = paths[ "summ" ][ "log_file" ]
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


def undistort( paths, conf ):
	"""Uses the fieldmaps to unwarp the functional images and create a mean image
	of all the unwarped functional images.
	"""

	func_fmap = paths[ "fmap" ][ "fmap_files" ][ 0 ]

	# motion-corrected images (input)
	func_corr = paths[ "func" ][ "corr_files" ]
	# unwarped images (output)
	func_uw = paths[ "func" ][ "uw_files" ]

	for i_run in xrange( len( func_corr )  ):

		fmri_tools.preproc.unwarp( func_corr[ i_run ],
		                           func_fmap,
		                           func_uw[ i_run ],
		                           conf[ "acq" ][ "dwell_ms" ],
		                           conf[ "acq" ][ "ph_encode_dir" ],
		                           log_path = paths[ "summ" ][ "log_file" ],
		                           pass_nocheck = False
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



def sess_reg( paths, conf ):
	"""Coregisters a within-session anatomical"""

	# need to pass the algorithm three translation parameters to get it close
	# to the mean functional

	par = []

	shft_max = 5

	for ( i_nudge, nudge_val ) in enumerate( conf[ "subj" ][ "nudge_vals" ] ):

		par.extend( [ "-parini",
		              "%d" % ( i_nudge + 1 ),
		              "%.3f" % nudge_val
		            ]
		          )

		par.extend( [ "-parang",
		              "%d" % ( i_nudge + 1 ),
		              "%.3f" % ( nudge_val - shft_max ),
		              "%.3f" % ( nudge_val + shft_max )
		            ]
		          )

	fmri_tools.preproc.anat_reg( paths[ "reg" ][ "base_dir" ],
	                             paths[ "reg" ][ "anat_file" ],
	                             paths[ "reg" ][ "mean_file" ],
	                             log_path = paths[ "summ" ][ "log_file" ],
	                             extra_allineate_opts = " ".join( par ),
	                             max_shft = None,
	                             max_rot = 10
	                           )


def vol_to_surf( paths, _ ):
	"""Converts the functional volume-based images to SUMA surfaces."""

	# this puts some output in the working directory, so change to where we want
	# to save stuff
	start_dir = os.getcwd()

	# images to project (unwarped)
	vol_files = paths[ "func" ][ "uw_files" ]
	# surface files to write
	surf_files = paths[ "func" ][ "surf_files" ]

	# number of increments to divide the surface interval
	f_steps = 15
	# how to deal with multiple nodes lying on a given voxel
	f_index = "nodes"

	for ( vol_file, surf_file ) in zip( vol_files, surf_files ):

		( file_dir, _ ) = os.path.split( vol_file )

		os.chdir( file_dir )

		for hemi in [ "lh", "rh" ]:

			surf_cmd = [ "3dVol2Surf",
			             "-spec", "%s%s.spec" % ( paths[ "reg" ][ "spec_base" ],
			                                      hemi
			                                    ),
			             "-surf_A", "smoothwm",
			             "-surf_B", "pial",
			             "-map_func", "ave",
			             "-f_steps", "%d" % f_steps,
			             "-f_index", f_index,
			             "-sv", "%s+orig" % paths[ "reg" ][ "reg_file" ],
			             "-grid_parent", "%s.nii" % vol_file,
			             "-out_niml", "%s_%s.niml.dset" % ( surf_file, hemi ),
			             "-overwrite"
			           ]

			fmri_tools.utils.run_cmd( surf_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

	os.chdir( start_dir )
