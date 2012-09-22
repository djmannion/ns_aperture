"""
Set of routines to pre-process the fMRI data for the natural scenes aperture
fMRI experiment.
"""

from __future__ import division

import os.path

import fmri_tools.preproc, fmri_tools.utils


def convert( paths, conf ):
	"""Converts the functionals and fieldmaps from dicom to nifti"""

	# functional DICOM files
	dcm_files = [ os.path.join( dcm_dir, dcm_file )
	              for ( dcm_dir, dcm_file ) in zip( paths[ "func" ][ "raw_dirs" ],
	                                                conf[ "subj" ][ "func_dcm" ]
	                                              )
	            ]

	# add in the fieldmap DICOMs
	dcm_files += [ os.path.join( paths[ "fmap" ][ "base_dir" ],
	                             "f01",
	                             "%s-raw" % im_type,
	                             dcm_file
	                           )
	               for ( im_type, dcm_file ) in zip( [ "mag", "ph" ],
	                                                 conf[ "subj" ][ "fmap_dcm" ]
	                                               )
	             ]

	# aggregate the images paths
	nii_files = ( paths[ "func" ][ "orig_files" ] +
	              paths[ "fmap" ][ "mag_files" ] +
	              paths[ "fmap" ][ "ph_files" ]
	            )

	for ( dcm_file, nii_file ) in zip( dcm_files, nii_files ):

		fmri_tools.preproc.dcm_to_nii( dcm_file,
		                               "%s.nii" % nii_file,
		                               reorient_dim = conf[ "acq" ][ "ras" ],
		                               log_path = paths[ "summ" ][ "log_file" ]
		                             )

	# generate the full paths (with assumed extension) of the newly-created nifti
	# files
	full_nii_paths = [ "%s.nii" % nii_file for nii_file in nii_files ]

	# check that they are all unique
	assert( fmri_tools.utils.files_are_unique( full_nii_paths ) )

	# files to go into the summary
	summ_paths = paths[ "func" ][ "orig_files" ]

	# make a summary image from the files
	fmri_tools.preproc.gen_sess_summ_img( summ_paths,
	                                      paths[ "summ" ][ "orig_summ_file" ],
	                                      log_path = paths[ "summ" ][ "log_file" ]
	                                    )


def motion_correct( paths, conf ):
	"""Performs motion correction"""

	# the index to the run that will be the 'base', to be corrected to
	# it is stored in one-based, hence the minus 1 to get to an index
	i_mc_base = conf[ "subj" ][ "mot_base" ] - 1

	mc_base = "%s.nii[0]" % paths[ "func" ][ "orig_files" ][ i_mc_base ]

	for ( orig_file, corr_file ) in zip( paths[ "func" ][ "orig_files" ],
	                                     paths[ "func" ][ "corr_files" ]
	                                   ):

		mc_cmd = [ "3dvolreg",
		           "-twopass",
		           "-prefix", "%s.nii" % corr_file,
		           "-overwrite",
		           "-base", mc_base,
		           "-zpad", "5",
		           "-heptic",  # fourier can cause ringing artefacts
		           "%s.nii" % orig_file
		         ]

		fmri_tools.utils.run_cmd( mc_cmd,
		                          env = fmri_tools.utils.get_env(),
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



def surf_reg( paths, conf ):
	"""Coregisters an anatomical with the SUMA reference"""

	base_anat = paths[ "reg" ][ "anat" ]
	reg_anat = paths[ "reg" ][ "reg_anat" ]

	# use the mean to represent the functional images
	base_func = "%s.nii" % paths[ "summ" ][ "mean_file" ]

	coreg_cmd = [ "3dAllineate",
	              "-base", base_func,
	              "-source", base_anat,
	              "-prefix", reg_anat,
	              "-cost", "nmi",  # normalised mutual info cost function
	              "-master", "SOURCE",
	              "-maxrot", "15",
	              "-overwrite",
	              "-warp", "shift_rotate",  # only rigid transforms
	              "-onepass",  # false minima if it is allowed to wander
	              "-verb"
	            ]

	# pass the algorithm three translation parameters to get it close to the
	# mean functional
	for ( i_nudge, nudge_val ) in enumerate( conf[ "subj" ][ "nudge_vals" ] ):

		coreg_cmd.extend( [ "-parini",
		                    "%d" % ( i_nudge + 1 ),
		                    "%.3f" % nudge_val
		                  ]
		                )

	fmri_tools.utils.run_cmd( coreg_cmd,
	                          env = fmri_tools.utils.get_env(),
	                          log_path = paths[ "summ" ][ "log_file" ]
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
			             "-spec", "%s%s.spec" % ( paths[ "reg" ][ "spec" ], hemi ),
			             "-surf_A", "smoothwm",
			             "-surf_B", "pial",
			             "-map_func", "ave",
			             "-f_steps", "%d" % f_steps,
			             "-f_index", f_index,
			             "-sv", paths[ "reg" ][ "reg_anat" ],
			             "-grid_parent", "%s.nii" % vol_file,
			             "-out_niml", "%s_%s.niml.dset" % ( surf_file, hemi ),
			             "-overwrite"
			           ]

			fmri_tools.utils.run_cmd( surf_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

	os.chdir( start_dir )
