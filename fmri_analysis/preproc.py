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


def exp_design_prep( paths, conf ):
	"""Prepares the designs for GLM analysis"""

	n_vols = int( conf[ "exp" ][ "run_len_s" ] / conf[ "acq" ][ "tr_s" ] )

	seq_info = ns_aperture.fmri.exp.get_seq_ind()

	# coherent blocks, initial coherent block, initial non-coherent block
	n_cond = 3

	cond_files = [ open( "%s%d.txt" % ( paths[ "ana" ][ "time_files" ],
	                                    cond_num
	                                  ),
	                     "w"
	                   )
	               for cond_num in np.arange( 1, n_cond + 1 )
	             ]

	for i_run in xrange( len( conf[ "subj" ][ "exp_runs" ] ) ):

		run_times = []
		run_conds = []

		run_seq = np.load( "%s%d.npy" % ( paths[ "log" ][ "seq_base" ],
		                                  i_run + 1
		                                )
		                 )

		( n_evt, n_params ) = run_seq.shape

		for i_evt in xrange( n_evt ):

			curr_block_num = run_seq[ i_evt, seq_info[ "block_num" ] ]
			prev_block_num = run_seq[ i_evt - 1, seq_info[ "block_num" ] ]

			is_transition = ( curr_block_num != prev_block_num )

			if is_transition:

				start_time_s = run_seq[ i_evt, seq_info[ "time_s" ] ]

				cond = int( run_seq[ i_evt, seq_info[ "block_type" ] ] )

				# don't want to model non-coherent blocks (for non-first blocks,
				# anyway)
				if ( cond == 1 ) and ( curr_block_num > 1 ):
					cond = 999

				if curr_block_num == 1:
					# this makes the first starting coherent block 1, first starting
					# non-coherent block 2
					cond += 1

				if cond != 999:
					run_times.append( start_time_s )
					run_conds.append( cond )

		run_times = np.array( run_times )
		run_conds = np.array( run_conds )

		for i_cond in xrange( n_cond ):

			i_evt_cond = np.where( run_conds == i_cond )[ 0 ]

			if i_evt_cond.size == 0:
				cond_files[ i_cond ].write( "*" )
			else:
				_ = [ cond_files[ i_cond ].write( "%.5f\t" % evt_time )
				      for evt_time in run_times[ i_evt_cond ]
				    ]

			cond_files[ i_cond ].write( "\n" )

	_ = [ cond_file.close() for cond_file in cond_files ]

	# POLYNOMIALS
	# ---

	n_pre_vol = int( conf[ "ana" ][ "exp_pre_cull_s" ] /
	                 conf[ "acq" ][ "tr_s" ]
	               )

	n_post_vol = int( conf[ "ana" ][ "exp_pre_cull_s" ] /
	                  conf[ "acq" ][ "tr_s" ]
	                )

	n_valid_vol = n_vols - n_pre_vol - n_post_vol

	# compute the polynomial timecourses
	run_trends = fmri_tools.utils.legendre_poly( conf[ "ana" ][ "poly_ord" ],
	                                             int( n_valid_vol ),
	                                             pre_n = n_pre_vol,
	                                             post_n = n_post_vol
	                                           )

	assert( run_trends.shape[ 0 ] == n_vols )

	n_runs = len( conf[ "subj" ][ "exp_runs" ] )

	# need to have a set of trends for each run, zeroed elsewhere
	bl_trends = np.zeros( ( n_vols * n_runs,
	                        conf[ "ana" ][ "poly_ord" ] * n_runs
	                      )
	                    )

	for i_run in xrange( n_runs ):

		i_row_start = i_run * run_trends.shape[ 0 ]
		i_row_end = i_row_start + run_trends.shape[ 0 ]

		i_col_start = i_run * run_trends.shape[ 1 ]
		i_col_end = i_col_start + run_trends.shape[ 1 ]

		bl_trends[ i_row_start:i_row_end, i_col_start:i_col_end ] = run_trends

	np.savetxt( paths[ "ana" ][ "bl_poly" ], bl_trends )


	# MOTION PARAMETERS
	# ---

	all_mc = np.loadtxt( paths[ "summ" ][ "mot_est_file" ] )

	exp_mc = all_mc[ :( n_vols * n_runs ), : ]

	np.savetxt( paths[ "ana" ][ "mot_est" ], exp_mc )


def loc_design_prep( paths, conf ):
	"""Prepares the designs for GLM analysis"""

	n_vols = int( conf[ "exp" ][ "loc_run_full_len_s" ] /
	              conf[ "acq" ][ "tr_s" ]
	            )

	seq_info = ns_aperture.fmri.loc.get_seq_ind()

	# L /R
	n_cond = 2

	cond_files = [ open( "%s%d.txt" % ( paths[ "loc" ][ "time_files" ],
	                                    cond_num
	                                  ),
	                     "w"
	                   )
	               for cond_num in np.arange( 1, n_cond + 1 )
	             ]

	loc_ord = ( "AB", "BA" )

	for i_run in xrange( conf[ "subj" ][ "n_loc_runs" ] ):

		run_times = []
		run_conds = []

		run_seq = ns_aperture.fmri.loc.get_seq( conf, loc_ord[ i_run ] )

		( n_evt, n_params ) = run_seq.shape

		for i_evt in xrange( n_evt ):

			curr_block_num = run_seq[ i_evt, seq_info[ "block_num" ] ]
			prev_block_num = run_seq[ i_evt - 1, seq_info[ "block_num" ] ]

			is_transition = ( curr_block_num != prev_block_num )

			if is_transition:

				start_time_s = run_seq[ i_evt, seq_info[ "time_s" ] ]

				cond = int( run_seq[ i_evt, seq_info[ "block_type" ] ] )

				# 0 = blank
				if cond > 0:
					run_times.append( start_time_s )
					run_conds.append( cond - 1 )

		run_times = np.array( run_times )
		run_conds = np.array( run_conds )

		for i_cond in xrange( n_cond ):

			i_evt_cond = np.where( run_conds == i_cond )[ 0 ]

			if i_evt_cond.size == 0:
				cond_files[ i_cond ].write( "*" )
			else:
				_ = [ cond_files[ i_cond ].write( "%.5f\t" % evt_time )
				      for evt_time in run_times[ i_evt_cond ]
				    ]

			cond_files[ i_cond ].write( "\n" )

	_ = [ cond_file.close() for cond_file in cond_files ]


	# MOTION PARAMETERS
	# ---

	all_mc = np.loadtxt( paths[ "summ" ][ "mot_est_file" ] )

	i_loc_start = ( conf[ "exp" ][ "run_len_s" ] /
	                conf[ "acq" ][ "tr_s" ] *
	                conf[ "subj" ][ "n_exp_runs" ]
	              )

	loc_mc = all_mc[ i_loc_start:, : ]

	assert( loc_mc.shape[ 0 ] == ( n_vols * conf[ "subj" ][ "n_loc_runs" ] ) )

	np.savetxt( paths[ "loc" ][ "mot_est" ], loc_mc )


def filter_for_svm( paths, conf ):
	"""Filters the timecourses for SVM analysis"""

	start_dir = os.getcwd()

	os.chdir( paths[ "svm" ][ "base_dir" ] )

	i_surf_files = np.array( conf[ "subj" ][ "exp_runs" ] ).astype( "int" ) - 1

	surf_files = [ paths[ "func" ][ "surf_files" ][ i_surf ]
	               for i_surf in i_surf_files
	             ]

	filt_files = [ paths[ "svm" ][ "filt_files" ][ i_surf ]
	               for i_surf in i_surf_files
	             ]

	for hemi in [ "lh", "rh" ]:

		for ( surf_file, filt_file ) in zip( surf_files, filt_files ):

			filt_cmd = [ "3dDetrend",
			             "-prefix", "%s_%s.niml.dset" % ( filt_file, hemi ),
			             "-polort", "%d" % ( conf[ "ana" ][ "poly_ord" ] - 1 ),
			             "-overwrite",
			             "%s_%s.niml.dset" % ( surf_file, hemi )
			           ]

			fmri_tools.utils.run_cmd( filt_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

			# convert to full
			full_filt_file = "%s_%s-full.niml.dset" % ( filt_file, hemi )
			pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

			fmri_tools.utils.sparse_to_full( "%s_%s.niml.dset" % ( filt_file, hemi ),
			                                 full_filt_file,
			                                 pad_node = pad_node,
			                                 log_path = paths[ "summ" ][ "log_file" ],
			                                 overwrite = True
			                               )

	os.chdir( start_dir )
