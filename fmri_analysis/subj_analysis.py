"""
Set of routines to analyse single-subject fMRI data for the natural scenes
apertures fMRI experiment.
"""

from __future__ import division

import os, os.path
import tempfile

import numpy as np
import scipy.stats

import fmri_tools.utils

import ns_aperture.fmri.exp

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

	# minus one because the range is inclusive
	pre_cens = conf[ "ana" ][ "exp_pre_cull_s" ] / conf[ "acq" ][ "tr_s" ] - 1
	post_cens = ( ( conf[ "exp" ][ "run_len_s" ] -
	                conf[ "ana" ][ "exp_post_cull_s" ]
	              ) /
	              conf[ "acq" ][ "tr_s" ] -
	              1
	            )

	# in AFNI-aware format; (runs):start-end
	censor_str = "*:0-%d,*:%d-%d" % ( pre_cens,
	                                  post_cens,
	                                  conf[ "exp" ][ "run_len_s" ] /
	                                  conf[ "acq" ][ "tr_s" ] - 1
	                                )

	for hemi in [ "lh", "rh" ]:

		glm_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "glm" ], hemi )
		beta_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "beta" ], hemi )

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
		                  "-ortvec", paths[ "ana" ][ "mot_est" ], "mot",
		                  "-local_times",
		                  "-CENSORTR", censor_str,
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-cbucket", beta_file,
		                  "-bucket", glm_file,
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


	os.chdir( start_dir )


def loc_glm( paths, conf ):
	"""Localiser GLM"""

	start_dir = os.getcwd()

	os.chdir( paths[ "loc" ][ "base_dir" ] )

	n_cond = 2

	hrf_model = conf[ "ana" ][ "hrf_model" ]

	stim_files = [ "%s%d.txt" % ( paths[ "loc" ][ "time_files" ],
	                              cond_num
	                            )
	               for cond_num in np.arange( 1, n_cond + 1 )
	             ]

	stim_labels = ( "L", "R" )

	con_str = [ "SYM: +L \ +R",  # F (OR) contrast
	            "SYM: +L -R"     # t (AND) contrast
	          ]

	con_lbl = [ "LorR", "LgtR" ]

	for hemi in [ "lh", "rh" ]:

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		i_surf_files = np.array( conf[ "subj" ][ "loc_runs" ] ).astype( "int" ) - 1

		surf_files = [ paths[ "func" ][ "surf_files" ][ i_surf ]
		               for i_surf in i_surf_files
		             ]

		glm_cmd.extend( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                  for surf_file in surf_files
		                ]
		              )

		glm_cmd.extend( [ "-force_TR", "%.3f" % conf[ "acq" ][ "tr_s" ],
		                  "-polort", "%d" % ( conf[ "ana" ][ "poly_ord" ] - 1 ),
		                  "-ortvec", paths[ "loc" ][ "mot_est" ], "mot",
		                  "-local_times",
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-overwrite",
		                  "-x1D_stop",  # want to use REML, so don't bother running
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

		# loop through all our contrasts
		for i_con in xrange( len( con_str ) ):

			glm_cmd.extend( [ "-gltsym",
			                  con_str[ i_con ]
			                ]
			              )

			glm_cmd.extend( [ "-glt_label",
			                  "%d" % ( i_con + 1 ),
			                  con_lbl[ i_con ]
			                ]
			              )

		# run this first GLM
		fmri_tools.utils.run_cmd( glm_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# delete the annoying command file that 3dDeconvolve writes
		os.remove( "Decon.REML_cmd" )

		glm_file = "%s_%s" % ( paths[ "loc" ][ "glm" ], hemi )
		beta_file = "%s_%s" % ( paths[ "loc" ][ "beta" ], hemi )

		reml_cmd = [ "3dREMLfit",
		             "-matrix", "exp_design.xmat.1D",
		             "-Rbeta", "%s_reml.niml.dset" % beta_file,
		             "-tout",
		             "-fout",
		             "-Rbuck", "%s_reml.niml.dset" % glm_file,
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

		# convert to full
		full_glm_file = "%s_%s_reml-full.niml.dset" % ( paths[ "loc" ][ "glm" ], hemi )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		fmri_tools.utils.sparse_to_full( "%s_reml.niml.dset" % glm_file,
		                                 full_glm_file,
		                                 pad_node = pad_node,
		                                 log_path = paths[ "summ" ][ "log_file" ],
		                                 overwrite = True
		                               )

	os.chdir( start_dir )



def loc_mask( paths, conf ):
	"""Creates a mask from the localiser contrast"""

	start_dir = os.getcwd()

	os.chdir( paths[ "loc" ][ "base_dir" ] )

	# bricks in the glm file that contains the localiser statistics
	# [ rvf T, lvf T ]
	loc_stat_bricks = [ 5, 2 ]

	for ( i_hemi, hemi ) in enumerate( [ "lh", "rh" ] ):

		glm_file = "%s_%s_reml.niml.dset" % ( paths[ "loc" ][ "glm" ], hemi )

		fdr_file = "%s_%s.niml.dset" % ( paths[ "loc" ][ "fdr" ], hemi )

		# convert the statistics for the localiser to a q (FDR) value
		fdr_cmd = [ "3dFDR",
		            "-input", glm_file,
		            "-prefix", fdr_file,
		            "-qval",  # specify that we want q, not z
		            "-float",
		            "-overwrite"
		          ]

		fmri_tools.utils.run_cmd( fdr_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		loc_mask_file = "%s_%s.niml.dset" % ( paths[ "loc" ][ "mask" ], hemi )

		q_thresh = conf[ "ana" ][ "q_thr" ]

		i_contra = i_hemi
		i_ipsi = np.logical_not( i_hemi )

		# 'a'
		contra_vf_p = "%s[%d]" % ( fdr_file, loc_stat_bricks[ i_contra ] )

		# 'b'
		ipsi_vf_p = "%s[%d]" % ( fdr_file, loc_stat_bricks[ i_ipsi ] )

		# 'c'
		contra_vf_t = "%s[%d]" % ( glm_file, loc_stat_bricks[ i_contra ] )

		expr = "and( within( a, 0, %(q).6f ), "  # contra p < thresh
		expr += "within( b, %(q).6f, 1 ), "  # ipsi p > thresh
		expr += "ispositive( c ) )"

		expr = expr % { "q" : q_thresh }

		# create a localiser mask as nodes that both have a q that is below
		# threshold and have positive beta weights
		mask_cmd = [ "3dcalc",
		             "-a", contra_vf_p,
		             "-b", ipsi_vf_p,
		             "-c", contra_vf_t,
		             "-expr", expr,
		             "-prefix", loc_mask_file,
		             "-overwrite"
		           ]

		fmri_tools.utils.run_cmd( mask_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# convert to full
		full_mask_file = "%s_%s-full" % ( paths[ "loc" ][ "mask" ], hemi )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		fmri_tools.utils.sparse_to_full( loc_mask_file,
		                                 full_mask_file,
		                                 pad_node = pad_node,
		                                 log_path = paths[ "summ" ][ "log_file" ],
		                                 overwrite = True
		                               )

	os.chdir( start_dir )


def beta_to_psc( paths, conf ):
	"""Convert the GLM beta weights into units of percent signal change"""

	# these are the indices into the beta files for the data we want to convert
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

	for hemi in [ "lh", "rh" ]:

		# the *full* localiser mask file
		loc_mask_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "loc_mask" ],
		                                           hemi
		                                         )

		# expression to apply the localiser mask
		cmask_expr = "-a %s -expr step(a)" % loc_mask_file

		# the *full* ROI file
		roi_file = "%s_%s-full.niml.dset" % ( paths[ "rois" ][ "dset" ], hemi )

		# iterate over all the ROIs
		for ( roi_name, roi_val ) in conf[ "ana" ][ "rois" ]:

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
			            "-mask", roi_file,
			            "-cmask", cmask_expr,
			            "-mrange", roi_val, roi_val,
			            "-noijk",
			            "-o", roi_psc_file,
			            data_file
			          ]

			fmri_tools.utils.run_cmd( xtr_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )


def raw_adj( paths, conf ):
	"""Concatenates raw timecourses and adjusts them for baselines"""

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "base_dir" ] )

	# since we use censoring in the GLM rather than trimming the data, we need to
	# specify a volume range here
	censor_vols = conf[ "exp" ][ "pre_len_s" ] / conf[ "acq" ][ "tr_s" ]
	brick_range = "[%d..$]" % censor_vols

	for hemi in [ "lh", "rh" ]:

		# create a raw input dataset, concatentated across all runs
		surf_files = [ "%s_%s.niml.dset%s" % ( surf_file, hemi, brick_range )
		               for surf_file in paths[ "func" ][ "surf_files" ]
		             ]

		raw_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "raw" ], hemi )

		cat_cmd = [ "3dTcat",
		            "-overwrite",
		            "-prefix", raw_file
		          ]

		cat_cmd.extend( surf_files )

		fmri_tools.utils.run_cmd( cat_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# adjust the raw file by subtracting the baseline
		raw_adj_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "raw_adj" ], hemi )

		bltc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "bltc" ], hemi )

		adj_cmd = [ "3dcalc",
		            "-fscale",
		            "-a", bltc_file,
		            "-b", raw_file,
		            "-expr", "b-a",
		            "-prefix", raw_adj_file,
		            "-overwrite"
		          ]

		fmri_tools.utils.run_cmd( adj_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# dataset holding the beta weights
		beta_file = "%s_%s_reml.niml.dset" % ( paths[ "ana" ][ "beta" ], hemi )

		# design matrix file
		mat_file = os.path.join( paths[ "ana" ][ "base_dir" ],
		                         "exp_design.xmat.1D"
		                       )

		pred_adj_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "pred_adj" ], hemi )

		# generate an signal timecourse
		bl_cmd = [ "3dSynthesize",
		           "-cbucket", beta_file,
		           "-matrix", mat_file,
		           "-cenfill", "none",  # important to match the censored data
		           "-select", "allfunc",
		           "-prefix", pred_adj_file,
		           "-overwrite"
		         ]

		fmri_tools.utils.run_cmd( bl_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# convert to full
		full_raw_file = "%s_%s-full" % ( paths[ "ana" ][ "raw" ], hemi )
		full_raw_adj_file = "%s_%s-full" % ( paths[ "ana" ][ "raw_adj" ], hemi )
		full_pred_adj_file = "%s_%s-full" % ( paths[ "ana" ][ "pred_adj" ], hemi )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		pad_files = [ [ raw_file, full_raw_file ],
		              [ raw_adj_file, full_raw_adj_file ],
		              [ pred_adj_file, full_pred_adj_file ]
		            ]

		for ( sparse_file, full_file ) in pad_files:

			fmri_tools.utils.sparse_to_full( sparse_file,
			                                 full_file,
			                                 pad_node = pad_node,
			                                 log_path = paths[ "summ" ][ "log_file" ],
			                                 overwrite = True
			                               )

	os.chdir( start_dir )


def svm_roi_xtr( paths, conf ):
	"""Extract the timecourses for each node of each ROI"""

	for ( roi_name, roi_val ) in conf[ "ana" ][ "rois" ]:

		roi_data = []

		for filt_file in paths[ "svm" ][ "filt_files" ]:

			run_data = []

			for hemi in [ "lh", "rh" ]:

				# the *full* ROI file
				roi_file = "%s_%s-full.niml.dset" % ( paths[ "rois" ][ "dset" ], hemi )

				out_file = os.path.join( paths[ "svm" ][ "base_dir" ], "temp.txt" )

				# our input dataset
				data_file = "%s_%s-full.niml.dset" % ( filt_file, hemi )

				# use the ROI file to mask the input dataset
				xtr_cmd = [ "3dmaskdump",
				            "-mask", roi_file,
				            "-mrange", roi_val, roi_val,
				            "-noijk",
				            "-o", out_file,
				            data_file
				          ]

				fmri_tools.utils.run_cmd( xtr_cmd,
				                          env = fmri_tools.utils.get_env(),
				                          log_path = paths[ "summ" ][ "log_file" ]
				                        )

				out_data = np.loadtxt( out_file )

				run_data.append( out_data )

				os.remove( out_file )

			run_data = np.concatenate( run_data ).T

			roi_data.append( run_data )

		# runs x time x node
		roi_data = np.array( roi_data )

		roi_data_file = "%s-%s.npy" % ( paths[ "svm" ][ "orig" ], roi_name )

		np.save( roi_data_file, roi_data )


def svm_roi_loc_stat( paths, conf ):
	"""Extract the timecourses for each node of each ROI"""

	start_dir = os.getcwd()

	os.chdir( paths[ "loc" ][ "base_dir" ] )

	# lh, rh
	stat_bricks = [ 2, 5 ]

	for ( roi_name, roi_val ) in conf[ "ana" ][ "rois" ]:

		hemi_data = []

		for hemi in [ "lh", "rh" ]:

			loc_file = "%s_%s_reml-full.niml.dset" % ( paths[ "loc" ][ "glm" ], hemi )
			comb_file = "%s_%s-full.niml.dset" % ( paths[ "loc" ][ "stat" ], hemi )

			# take the maximum of the two conjuction contrasts
			stat_cmd = [ "3dcalc",
			             "-a", "%s[%d]" % ( loc_file, stat_bricks[ 0 ] ),
			             "-b", "%s[%d]" % ( loc_file, stat_bricks[ 1 ] ),
			             "-expr", "max( a, b )",
			             "-prefix", comb_file,
			             "-overwrite"
			           ]

			fmri_tools.utils.run_cmd( stat_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

			# the *full* ROI file
			roi_file = "%s_%s-full.niml.dset" % ( paths[ "rois" ][ "dset" ], hemi )

			out_file = os.path.join( paths[ "svm" ][ "base_dir" ], "temp.txt" )

			xtr_cmd = [ "3dmaskdump",
			            "-mask", roi_file,
			            "-mrange", roi_val, roi_val,
			            "-noijk",
			            "-o", out_file,
			            comb_file
			          ]

			fmri_tools.utils.run_cmd( xtr_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

			out_data = np.loadtxt( out_file )

			hemi_data.append( out_data )

			os.remove( out_file )

		hemi_data = np.concatenate( hemi_data ).T

		#  x node
		roi_data = np.array( hemi_data )

		roi_data_file = "%s-%s.npy" % ( paths[ "svm" ][ "loc_stat" ], roi_name )

		np.save( roi_data_file, roi_data )

	os.chdir( start_dir )


def svm_info( paths, conf ):
	"""Extract z-scored blocks for each ROI"""

	seq_info = ns_aperture.fmri.exp.get_seq_ind()

	run_info = np.empty( ( conf[ "subj" ][ "n_exp_runs" ],
	                       conf[ "exp" ][ "n_blocks" ] - 2,
	                       2  # onset volume, block type
	                     )
	                   )
	run_info.fill( np.NAN )

	for i_run in xrange( len( conf[ "subj" ][ "exp_runs" ] ) ):

		run_seq = np.load( "%s%d.npy" % ( paths[ "log" ][ "seq_base" ],
		                                  i_run + 1
		                                )
		                 )

		for i_evt in xrange( run_seq.shape[ 0 ] ):

			curr_block_num = run_seq[ i_evt, seq_info[ "block_num" ] ]
			prev_block_num = run_seq[ i_evt - 1, seq_info[ "block_num" ] ]

			is_transition = ( curr_block_num != prev_block_num )

			if is_transition:

				if curr_block_num != 1 and curr_block_num != conf[ "exp" ][ "n_blocks" ]:

					cond = run_seq[ i_evt, seq_info[ "block_type" ] ]

					if cond == 0:
						block_type = 1
					else:
						block_type = -1

					i_block = curr_block_num - 1 - 1

					run_info[ i_run, i_block, 1 ] = block_type

					time_s = run_seq[ i_evt, seq_info[ "time_s" ] ]

					time_vol = time_s / conf[ "acq" ][ "tr_s" ]

					time_adj = time_vol + conf[ "ana" ][ "hrf_corr_vol" ]

					run_info[ i_run, i_block, 0 ] = time_adj

	np.save( paths[ "svm" ][ "run_info" ], run_info )


def svm_roi_z( paths, conf ):
	"""a"""

	run_info = np.load( paths[ "svm" ][ "run_info" ] )

	n_vol_per_block = conf[ "exp" ][ "block_len_s" ] / conf[ "acq" ][ "tr_s" ]

	for ( roi_name, _ ) in conf[ "ana" ][ "rois" ]:

		# runs x time x nodes
		orig_data = np.load( "%s-%s.npy" % ( paths[ "svm" ][ "orig" ], roi_name ) )

		blk_data = np.empty( ( orig_data.shape[ 0 ],
		                       run_info.shape[ 1 ],
		                       orig_data.shape[ 2 ]
		                     )
		                   )
		blk_data.fill( np.NAN )

		for i_run in xrange( run_info.shape[ 0 ] ):

			for i_block in xrange( run_info.shape[ 1 ] ):

				i_r = np.arange( run_info[ i_run, i_block, 0 ],
				                 run_info[ i_run, i_block, 0 ] + n_vol_per_block
				               ).astype( "int" )

				# average over block timepoints
				data = np.mean( orig_data[ i_run, i_r, : ], axis = 0 )

				blk_data[ i_run, i_block, : ] = data

			blk_data[ i_run, :, : ] = scipy.stats.zscore( blk_data[ i_run, :, : ],
			                                              axis = 0
			                                            )

		roi_z_file = "%s-%s.npy" % ( paths[ "svm" ][ "z" ], roi_name )

		np.save( roi_z_file, blk_data )


def svm( paths, conf ):
	"""a"""

	train_runs = [ [ 2, 3, 4, 5, 6, 7, 8, 9 ],
	               [ 4, 5, 6, 7, 8, 9, 0, 1 ],
	               [ 6, 7, 8, 9, 0, 1, 2, 3 ],
	               [ 8, 9, 0, 1, 2, 3, 4, 5 ],
	               [ 0, 1, 2, 3, 4, 5, 6, 7 ]
	             ]

	test_runs = [ [ 0, 1 ],
	              [ 2, 3 ],
	              [ 4, 5 ],
	              [ 6, 7 ],
	              [ 8, 9 ]
	            ]

	svm_runs = zip( train_runs, test_runs )

	run_info = np.load( paths[ "svm" ][ "run_info" ] )

	node_range = np.arange( 10, 1001, 10 )

	for ( roi_name, _ ) in conf[ "ana" ][ "rois" ]:

		z = np.load( "%s-%s.npy" % ( paths[ "svm" ][ "z" ], roi_name ) )
		loc = np.load( "%s-%s.npy" % ( paths[ "svm" ][ "loc_stat" ], roi_name ) )

		i_loc_sort = np.argsort( loc )[ ::-1 ]

		svm_pred = np.empty( ( len( node_range ),  # nodes
		                       len( svm_runs ),    # folds
		                       len( test_runs[ 0 ] ) * z.shape[ 1 ],  # examples
		                       2  # true, predicted
		                     )
		                   )
		svm_pred.fill( np.NAN )

		for ( i_fold, ( train_set, test_set ) ) in enumerate( svm_runs ):

			fold_dir = os.path.join( "%s%02d" % ( paths[ "svm" ][ "fold_base" ],
			                                    i_fold + 1
			                                  )
			                       )

			for ( i_node_k, n_nodes ) in enumerate( node_range ):

				train_data = z[ train_set, :, : ]
				train_data = train_data[ :, :, i_loc_sort[ :n_nodes ] ]

				test_data = z[ test_set, :, : ]
				test_data = test_data[ :, :, i_loc_sort[ :n_nodes ] ]


				train_path = os.path.join( fold_dir,
				                           "%s_%s_fold_%02d_%04d_nodes.txt" % (
				                             paths[ "svm" ][ "train_base" ],
				                             roi_name,
				                             i_fold + 1,
				                             n_nodes
				                           )
				                         )

				train_file = open( train_path, "w" )

				for i_run in xrange( train_data.shape[ 0 ] ):
					for i_block in xrange( train_data.shape[ 1 ] ):

						ex_cond = run_info[ train_set[ i_run ], i_block, 1 ]

						ex_str = "%+d" % ex_cond

						for i_node in xrange( train_data.shape[ 2 ] ):

							ex_str += " %d:%.12f" % ( i_node + 1,
							                          train_data[ i_run, i_block, i_node ]
							                        )

						train_file.write( "%s\n" % ex_str )

				train_file.close()

				# do the training
				model_path = os.path.join( fold_dir,
				                           "%s_%s_fold_%02d_%04d_nodes.txt" % (
				                             paths[ "svm" ][ "model_base" ],
				                             roi_name,
				                             i_fold + 1,
				                             n_nodes
				                           )
				                         )

				train_cmd = [ "svm_learn",
				              train_path,
				              model_path
				            ]

				fmri_tools.utils.run_cmd( train_cmd,
				                          env = fmri_tools.utils.get_env(),
				                          log_path = paths[ "summ" ][ "log_file" ]
				                        )

				# write the test data file
				test_path = os.path.join( fold_dir,
				                          "%s_%s_fold_%02d_%04d_nodes.txt" % (
				                            paths[ "svm" ][ "test_base" ],
				                            roi_name,
				                            i_fold + 1,
				                            n_nodes
				                          )
				                        )

				test_file = open( test_path, "w" )

				i_count = 0

				for i_run in xrange( test_data.shape[ 0 ] ):
					for i_block in xrange( test_data.shape[ 1 ] ):

						ex_cond = run_info[ test_set[ i_run ], i_block, 1 ]

						svm_pred[ i_node_k, i_fold, i_count, 0 ] = ex_cond

						i_count += 1

						ex_str = "%+d" % ex_cond

						for i_node in xrange( test_data.shape[ 2 ] ):

							ex_str += " %d:%.12f" % ( i_node + 1,
							                          test_data[ i_run, i_block, i_node ]
							                        )

						test_file.write( "%s\n" % ex_str )

				test_file.close()

				# do the testing
				pred_path = os.path.join( fold_dir,
				                           "%s_%s_fold_%02d_%04d_nodes.txt" % (
				                             paths[ "svm" ][ "pred_base" ],
				                             roi_name,
				                             i_fold + 1,
				                             n_nodes
				                           )
				                         )

				test_cmd = [ "svm_classify",
				              test_path,
				              model_path,
				              pred_path
				            ]

				fmri_tools.utils.run_cmd( test_cmd,
				                          env = fmri_tools.utils.get_env(),
				                          log_path = paths[ "summ" ][ "log_file" ]
				                        )

				pred = np.loadtxt( pred_path )

				svm_pred[ i_node_k, i_fold, :, 1 ] = pred


		pred_summ_path = "%s_%s.npy" % ( paths[ "svm" ][ "summ" ], roi_name )

		np.save( pred_summ_path, svm_pred )


def roi_tc( paths, conf ):
	"""Compile raw and predicted (adjusted) timecourses for each ROI."""

	for hemi in [ "lh", "rh" ]:

		# the *full* localiser mask file
		loc_mask_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "loc_mask" ],
		                                           hemi
		                                         )

		# expression to apply the localiser mask
		cmask_expr = "-a %s -expr step(a)" % loc_mask_file

		# the *full* ROI file
		roi_file = "%s_%s-full.niml.dset" % ( paths[ "rois" ][ "dset" ], hemi )

		# iterate over all the ROIs
		for ( roi_name, roi_val ) in conf[ "ana" ][ "rois" ]:

			roi_raw_adj_file = "%s_%s_%s.txt" % ( paths[ "rois" ][ "raw_adj_tc" ],
			                                      roi_name,
			                                      hemi
			                                    )

			raw_adj_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "raw_adj" ],
			                                          hemi
			                                        )

			roi_pred_adj_file = "%s_%s_%s.txt" % ( paths[ "rois" ][ "pred_adj_tc" ],
			                                       roi_name,
			                                       hemi
			                                     )

			pred_adj_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "pred_adj" ],
			                                           hemi
			                                         )

			data_files = [ [ roi_raw_adj_file, raw_adj_file ],
			               [ roi_pred_adj_file, pred_adj_file ]
			             ]

			for ( out_file, in_file ) in data_files:

				# 3dmaskdump won't overwrite, so need to manually remove any prior data
				if os.path.exists( out_file ):
					os.remove( out_file )

				# use the ROI file to mask the input dataset
				xtr_cmd = [ "3dmaskdump",
				            "-mask", roi_file,
				            "-cmask", cmask_expr,
				            "-mrange", roi_val, roi_val,
				            "-noijk",
				            "-o", out_file,
				            in_file
				          ]

				fmri_tools.utils.run_cmd( xtr_cmd,
				                          env = fmri_tools.utils.get_env(),
				                          log_path = paths[ "summ" ][ "log_file" ]
				                        )
