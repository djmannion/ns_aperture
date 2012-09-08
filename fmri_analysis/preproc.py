"""
Set of routines to pre-process the fMRI data for the natural scenes aperture
fMRI experiment.
"""

from __future__ import division

import os.path
import tempfile

import numpy as np

import fmri_tools.preproc, fmri_tools.utils

import ns_aperture.fmri.exp, ns_aperture.fmri.loc


def convert( paths, conf ):
	"""Converts the functionals and fieldmaps from dicom to nifti"""

	# aggregate the dicom directories
	raw_dirs = ( paths[ "func" ][ "raw_dirs" ] +
	             paths[ "fmap" ][ "raw_mag_dirs" ] +
	             paths[ "fmap" ][ "raw_ph_dirs" ]
	           )

	# aggregate the output directories
	nii_dirs = ( paths[ "func" ][ "run_dirs" ] +
	             paths[ "fmap" ][ "fmap_dirs" ] +
	             paths[ "fmap" ][ "fmap_dirs" ]
	           )

	# aggregate the images paths
	img_paths = ( paths[ "func" ][ "orig_files" ] +
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

	mc_params = []

	for i_run in xrange( conf[ "subj" ][ "n_runs" ] ):

		orig_file = paths[ "func" ][ "orig_files" ][ i_run ]
		corr_file = paths[ "func" ][ "corr_files" ][ i_run ]

		# because we want to aggregate the motion correction files over the whole
		# session, we only store the individual run correction temporarily
		mc_txt = tempfile.NamedTemporaryFile()

		mc_cmd = [ "3dvolreg",
		           "-twopass",
		           "-prefix", "%s.nii" % corr_file,
		           "-1Dfile", mc_txt.name,
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

		# deal with the motion estimates
		mc_params.append( np.loadtxt( mc_txt.name ) )

	# concatenate the motion estimates over time
	mc_params = np.vstack( mc_params )

	np.savetxt( paths[ "summ" ][ "mot_est_file" ], mc_params )

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


def vol_to_surf( paths, conf ):
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


def design_prep( paths, conf ):
	"""Prepares the designs for GLM analysis"""

	# first, prepare the 'onset' regressor; this models the tail of the first
	# block response, which remains after excluding the timepoints corresponding
	# to the first block's stimulation
	exp_run_len_vol = conf[ "exp" ][ "run_len_s" ] / conf[ "acq" ][ "tr_s" ]

	x = tempfile.NamedTemporaryFile()

	onset_cmd = [ "3dDeconvolve",
	              "-polort", "-1",
	              "-nodata",  "%d" % exp_run_len_vol, "%.3f" % conf[ "acq" ][ "tr_s" ],
	              "-local_times",
	              "-num_stimts", "1",
	              "-stim_times", "1", "1D: 0", "SPMG1(16)",
	              "-x1D", x.name,
	              "-x1D_stop"
	            ]

	fmri_tools.utils.run_cmd( onset_cmd,
	                          env = fmri_tools.utils.get_env(),
	                          log_path = paths[ "summ" ][ "log_file" ]
	                        )

	tc = np.loadtxt( "%s.xmat.1D" % x.name )

	start_vol = int( conf[ "ana" ][ "exp_run_start_s" ] / conf[ "acq" ][ "tr_s" ] )
	n_vol = int( conf[ "ana" ][ "exp_run_dur_s" ] / conf[ "acq" ][ "tr_s" ] )

	tc = tc[ start_vol:( start_vol + n_vol ) ]

	A_reg = np.tile( np.hstack( ( tc, np.zeros( len( tc ) ) ) ),
	                 conf[ "exp" ][ "n_runs" ] / 2
	               )

	np.savetxt( paths[ "log" ][ "reg_A" ], A_reg, "%.16f" )

	B_reg = np.tile( np.hstack( ( np.zeros( len( tc ) ), tc ) ),
	                 conf[ "exp" ][ "n_runs" ] / 2
	               )

	np.savetxt( paths[ "log" ][ "reg_B" ], B_reg, "%.16f" )

	# exp
	run_file = open( paths[ "ana" ][ "exp_time_file" ], "w" )

	# get info on what index corresponds to what, in the run sequence file
	seq_ind = ns_aperture.fmri.exp.get_seq_ind()

	for i_run in xrange( conf[ "subj" ][ "n_runs" ] ):

		run_times = []

		# load the sequence for this run
		run_seq = np.load( "%s%d.npy" % ( paths[ "log" ][ "seq_base" ], i_run + 1 ) )

		# loop over each event in the sequence
		for i_evt in xrange( run_seq.shape[ 0 ] ):

			# pull out the event info
			evt_time_s = run_seq[ i_evt, seq_ind[ "time_s" ] ]
			evt_block_num = run_seq[ i_evt, seq_ind[ "block_num" ] ]
			evt_prev_block_num = run_seq[ i_evt - 1, seq_ind[ "block_num" ] ]
			evt_block_type = run_seq[ i_evt, seq_ind[ "block_type" ] ]
			evt_img_i_L = run_seq[ i_evt, seq_ind[ "img_i_L" ] ]
			evt_img_i_R = run_seq[ i_evt, seq_ind[ "img_i_R" ] ]

			# test if the event marks the first of a new block
			is_transition = ( evt_block_num != evt_prev_block_num )

			# test if the event / block is 'coherent'; if it's condition index is 0
			is_coh = ( evt_block_type == 0 )

			# only noteworthy if it is both the start of a new block and the block is
			# of coherent stimuli
			if np.logical_and( is_transition, is_coh ):

				# do a sanity check that it is in fact a coherent block by testing
				# whether both stimuli have the same id
				assert( evt_img_i_L == evt_img_i_R )

				run_times.append( evt_time_s )

		run_times = np.array( run_times )

		# subtract the time that we cull from the data
		run_times -= conf[ "ana" ][ "exp_run_start_s" ]

		# and remove any times that are outside of our data window
		ok = np.logical_and( run_times >= 0,
		                     run_times < conf[ "ana" ][ "exp_run_dur_s" ]
		                   )

		run_times = run_times[ ok ]

		# write the onset times to the run file
		for run_time in run_times:
			run_file.write( "%.5f\t" % run_time )

		run_file.write( "\n" )

	run_file.close()

	# loc
	loc_order = [ "AB", "BA" ]

	loc_time_files = [ open( loc_time_file, "w" )
	                   for loc_time_file in paths[ "ana" ][ "loc_time_files" ]
	                 ]

	for loc_ord in loc_order:

		block_seq = ns_aperture.fmri.loc.get_block_seq( loc_ord,
		                                                conf[ "exp" ][ "loc_n_blocks" ]
		                                              )

		i_lvf = np.where( block_seq == 1 )[ 0 ]
		lvf_t = i_lvf * conf[ "exp" ][ "block_len_s" ]

		for t in lvf_t:
			loc_time_files[ 0 ].write( "%.5f\t" % t )

		i_rvf = np.where( block_seq == 2 )[ 0 ]
		rvf_t = i_rvf * conf[ "exp" ][ "block_len_s" ]

		for t in rvf_t:
			loc_time_files[ 1 ].write( "%.5f\t" % t )

		loc_time_files[ 0 ].write( "\n" )
		loc_time_files[ 1 ].write( "\n" )

	loc_time_files[ 0 ].close()
	loc_time_files[ 1 ].close()

