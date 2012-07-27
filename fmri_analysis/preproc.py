"""
Set of routines to pre-process the fMRI data for the natural scenes aperture
fMRI experiment.
"""

from __future__ import division

import os.path
import tempfile

import numpy as np

import ns_aperture.fmri.exp, ns_aperture.fmri.loc
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
	                                      log_path = paths[ "summ" ][ "log_file" ]
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
	orig_paths = [ paths[ "func_%s" % im_type ][ "orig_files" ][ i_run - 1 ]
	               for i_run, im_type in run_order
	             ]
	corr_paths = [ paths[ "func_%s" % im_type ][ "corr_files" ][ i_run - 1 ]
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

	for i_run in xrange( len( func_corr )  ):

		fmri_tools.preproc.unwarp( func_corr[ i_run ],
		                           func_fmap[ i_run ],
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


def trim( paths, conf ):
	"""Trims the timecourses"""

	exp_start_vol = conf[ "ana" ][ "exp_run_start_s" ] / conf[ "acq" ][ "tr_s" ]
	exp_n_vol = conf[ "ana" ][ "exp_run_dur_s" ] / conf[ "acq" ][ "tr_s" ]

	for ( uw_file, trim_file ) in zip( paths[ "func_exp" ][ "uw_files" ],
	                                   paths[ "func_exp" ][ "trim_files" ]
	                                 ):

		exp_trim_cmd = [ "fslroi",
		                 uw_file,
		                 trim_file,
		                 "%d" % exp_start_vol,
		                 "%d" % exp_n_vol
		               ]

		fmri_tools.utils.run_cmd( exp_trim_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

	loc_start_vol = conf[ "ana" ][ "loc_run_start_s" ] / conf[ "acq" ][ "tr_s" ]
	loc_n_vol = conf[ "ana" ][ "loc_run_dur_s" ] / conf[ "acq" ][ "tr_s" ]

	for ( uw_file, trim_file ) in zip( paths[ "func_loc" ][ "uw_files" ],
	                                   paths[ "func_loc" ][ "trim_files" ]
	                                 ):

		loc_trim_cmd = [ "fslroi",
		                 uw_file,
		                 trim_file,
		                 "%d" % loc_start_vol,
		                 "%d" % loc_n_vol
		               ]

		fmri_tools.utils.run_cmd( loc_trim_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )


def surf_reg( paths, conf ):
	"""Coregisters an anatomical with the SUMA reference"""

	fmri_tools.preproc.surf_reg( paths[ "reg" ][ "rs_exp_anat" ],
	                             paths[ "reg" ][ "surf_anat" ],
	                             paths[ "summ" ][ "log_file" ]
	                           )


def vol_to_surf( paths, conf ):
	"""Converts the functional volume-based images to SUMA surfaces."""

	start_dir = os.getcwd()

	vol_files = ( paths[ "func_exp" ][ "trim_files" ] +
	              paths[ "func_loc" ][ "trim_files" ]
	            )

	surf_files = ( paths[ "func_exp" ][ "surf_files" ] +
	               paths[ "func_loc" ][ "surf_files" ]
	             )

	for ( vol_file, surf_file ) in zip( vol_files, surf_files ):

		file_dir = os.path.split( vol_file )[ 0 ]

		os.chdir( file_dir )

		for hemi in [ "lh", "rh" ]:

			out = "%s_%s.niml.dset" % ( surf_file,
			                            hemi
			                          )

			surf_cmd = [ "3dVol2Surf",
			             "-spec", "%s%s.spec" % ( paths[ "reg" ][ "spec" ], hemi ),
			             "-surf_A", "smoothwm",
			             "-surf_B", "pial",
			             "-map_func", "ave",
			             "-f_steps", "15",
			             "-f_index", "nodes",
			             "-sv", paths[ "reg" ][ "reg" ],
			             "-grid_parent", "%s.nii" % vol_file,
			             "-out_niml", out,
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

