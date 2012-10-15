"""
Set of routines to analyse single-subject fMRI data for the natural scenes
apertures fMRI experiment.
"""

from __future__ import division

import os, os.path

import numpy as np

import fmri_tools.utils

import ns_aperture.fmri.exp


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

		( n_evt, _ ) = run_seq.shape

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

	# CENSORING
	# ---

	cens = np.ones( n_vols )

	n_vols_per_block = int( conf[ "exp" ][ "block_len_s" ] /
	                        conf[ "acq" ][ "tr_s" ]
	                      )

	cens[ :n_vols_per_block ] = 0
	cens[ -n_vols_per_block: ] = 0

	assert( np.sum( cens == 0 ) == ( n_vols_per_block * 2 ) )

	cens = np.tile( cens, conf[ "subj" ][ "n_exp_runs" ] )

	np.savetxt( paths[ "ana" ][ "cens" ], cens )


def glm( paths, conf ):
	"""Experiment GLM"""

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "base_dir" ] )

	n_cond = 3

	hrf_model = conf[ "ana" ][ "hrf_model" ]

	stim_files = [ "%s%d.txt" % ( paths[ "ana" ][ "time_files" ],
	                              cond_num
	                            )
	               for cond_num in np.arange( 1, n_cond + 1 )
	             ]

	stim_labels = [ "coh", "Scoh", "Snoncoh" ]

	for hemi in [ "lh", "rh" ]:

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		i_surf_files = np.array( conf[ "subj" ][ "exp_runs" ] ).astype( "int" ) - 1

		surf_files = [ paths[ "func" ][ "surf_files" ][ i_surf ]
		               for i_surf in i_surf_files
		             ]

		glm_cmd.extend( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                  for surf_file in surf_files
		                ]
		              )

		glm_cmd.extend( [ "-force_TR", "%.3f" % conf[ "acq" ][ "tr_s" ],
		                  "-polort", "-1",  # we pass our own below
		                  "-ortvec", paths[ "ana" ][ "bl_poly" ], "poly",
		                  "-local_times",
		                  "-censor", paths[ "ana" ][ "cens" ],
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-x1D_stop",
		                  "-overwrite",
		                  "-num_stimts", "%d" % n_cond
		                ]
		              )

		for i_stim in xrange( n_cond ):

			glm_cmd.extend( [ "-stim_label",
			                  "%d" % ( i_stim + 1 ),
			                  stim_labels[ i_stim ]
			                ]
			              )

			glm_cmd.extend( [ "-stim_times",
			                  "%d" % ( i_stim + 1 ),
			                  stim_files[ i_stim ],
			                  hrf_model
			                ]
			              )

		# run the GLM
		fmri_tools.utils.run_cmd( glm_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# delete the annoying command file that 3dDeconvolve writes
		os.remove( "Decon.REML_cmd" )

		glm_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "glm" ], hemi )
		beta_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "beta" ], hemi )

		reml_cmd = [ "3dREMLfit",
		             "-matrix", "exp_design.xmat.1D",
		             "-tout",
		             "-Rbuck", glm_file,
		             "-Rbeta", beta_file,
		             "-overwrite",
		             "-input"
		           ]

		reml_cmd.append( " ".join( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                             for surf_file in surf_files
		                           ]
		                         )
		               )

		# run the proper GLM
		fmri_tools.utils.run_cmd( reml_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

	os.chdir( start_dir )


def beta_to_psc( paths, conf ):
	"""Convert the GLM beta weights into units of percent signal change"""

	# this is the index in the beta file for the data we want to convert
	beta_brick = "0"

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "base_dir" ] )

	for hemi in [ "lh", "rh" ]:

		# dataset holding the beta weights
		beta_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "beta" ], hemi )

		# design matrix file
		mat_file = os.path.join( paths[ "ana" ][ "base_dir" ],
		                         "exp_design.xmat.1D"
		                       )

		# baseline timecourse dataset, to write
		bltc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "bltc" ], hemi )

		# generate an average baseline timecourse
		bl_cmd = [ "3dSynthesize",
		           "-cbucket", beta_file,
		           "-matrix", mat_file,
		           "-cenfill", "none",  # this is important
		           "-select", "poly",  # only use the polynomials
		           "-prefix", bltc_file,
		           "-overwrite"
		         ]

		fmri_tools.utils.run_cmd( bl_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# baseline (point-estimate) dataset, to write
		bl_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "bl" ], hemi )

		# average baseline timecourse across time
		avg_cmd = [ "3dTstat",
		            "-mean",
		            "-overwrite",
		            "-prefix", bl_file,
		            bltc_file
		          ]

		fmri_tools.utils.run_cmd( avg_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# dataset to hold the percent signal change, to write
		psc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "psc" ], hemi )

		# the input beta file, with sub-brick selector
		beta_sel = "%s[%s]" % ( beta_file, beta_brick )

		# check that the label is as expected
		beta_label = fmri_tools.utils.get_dset_label( beta_sel )
		assert( beta_label == [ "coh#0" ] )

		# compute psc
		# from http://afni.nimh.nih.gov/sscc/gangc/TempNorm.html
		psc_cmd = [ "3dcalc",
		            "-fscale",
		            "-a", bl_file,
		            "-b", beta_sel,
		            "-expr", "100 * b/a * step (1- abs(b/a))",
		            "-prefix", psc_file,
		            "-overwrite"
		          ]

		fmri_tools.utils.run_cmd( psc_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# convert to full
		full_psc_file = "%s_%s-full" % ( paths[ "ana" ][ "psc" ], hemi )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		fmri_tools.utils.sparse_to_full( psc_file,
		                                 full_psc_file,
		                                 pad_node = pad_node,
		                                 log_path = paths[ "summ" ][ "log_file" ],
		                                 overwrite = True
		                               )

	os.chdir( start_dir )


def roi_xtr( paths, conf ):
	"""Extract PSC and statistics data from ROIs"""

	os.chdir( paths[ "rois" ][ "base_dir" ] )

	for ( roi_name, _ ) in conf[ "ana" ][ "rois" ]:

		for hemi in [ "lh", "rh" ]:

			# the *full* localiser mask file
			loc_mask_file = "%s_%s_%s-full.niml.dset" % ( paths[ "loc" ][ "roi_parc" ],
			                                              roi_name,
			                                              hemi
			                                            )

			roi_psc_file = "%s_%s_%s.txt" % ( paths[ "rois" ][ "psc" ],
			                                  roi_name,
			                                  hemi
			                                )

			# 3dmaskdump won't overwrite, so need to manually remove any prior file
			if os.path.exists( roi_psc_file ):
				os.remove( roi_psc_file )

			# our input dataset
			data_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "psc" ],
			                                       hemi
			                                     )

			# use the ROI file to mask the input dataset
			xtr_cmd = [ "3dmaskdump",
			            "-mask", loc_mask_file,
			            "-noijk",
			            "-o", roi_psc_file,
			            data_file,
			            loc_mask_file
			          ]

			fmri_tools.utils.run_cmd( xtr_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )


def roi_mean( paths, conf ):
	"""Average nodes within each parcel of each ROI"""

	os.chdir( paths[ "rois" ][ "base_dir" ] )

	for ( roi_name, _ ) in conf[ "ana" ][ "rois" ]:

		psc = np.vstack( [ np.loadtxt( "%s_%s_%s.txt" % ( paths[ "rois" ][ "psc" ],
		                                                  roi_name,
		                                                  hemi
		                                                )
		                             )
		                   for hemi in [ "lh", "rh" ]
		                 ]
		               )

		parc_mean = np.empty( len( conf[ "ana" ][ "i_parc" ] ) )
		parc_mean.fill( np.NAN )

		for ( i_parc, parc_ind ) in enumerate( conf[ "ana" ][ "i_parc" ] ):

			# psc array indices as those with a parcel number equal to any in the
			# list
			i_parc_psc = np.any( [ psc[ :, 1 ] == k
			                       for k in parc_ind
			                     ],
			                     axis = 0
			                   )

			parc_mean[ i_parc ] = np.mean( psc[ i_parc_psc, 0 ] )

		np.savetxt( "%s_%s.txt" % ( paths[ "rois" ][ "parc_psc" ],
		                            roi_name
		                          ),
		            parc_mean
		          )


def group_agg( grp_paths, conf ):
	"""Aggregate the univariate analyses across subjects"""

	n_subj = len( conf[ "all_subj" ] )

	n_rois = len( conf[ "ana" ][ "rois" ] ) * len( conf[ "ana" ][ "i_parc" ] )

	uni_data = np.empty( n_subj, n_rois )
	uni_data.fill( np.NAN )

	for ( i_subj, subj_id ) in conf[ "all_subj" ]:

		subj_conf = ns_aperture.config.get_conf( subj_id )
		subj_paths = ns_aperture.fmri_analysis.paths.get_subj_paths( subj_conf )

		i_roi = 0

		for ( roi_name, _ ) in conf[ "ana" ][ "rois" ]:

			roi_psc = np.loadtxt( paths[ "rois" ][ "parc_psc" ], roi_name )

			uni_data[ i_subj, i_roi:( i_roi + 2 ) ] = roi_psc

			i_roi += 2

	np.savetxt( grp_paths[ "uni" ], uni_data )
