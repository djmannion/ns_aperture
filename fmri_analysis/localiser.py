"""
Set of routines to analyse single-subject localiser fMRI data for the natural
scenes apertures fMRI experiment.
"""

from __future__ import division

import os, os.path
import tempfile

import numpy as np
import scipy.stats

import fmri_tools.utils

import ns_aperture.fmri.exp


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


def loc_parcellate( paths, conf ):
	"""Extract t statistics from rois"""

	# LVF, RVF
	loc_stat_bricks = "2,4"

	loc_rois = [ [ "v1", "1" ],
	             [ "v2", "2" ],
	             [ "v3", "3" ]
	           ]

	os.chdir( paths[ "loc" ][ "base_dir" ] )

	# check by running:
	#  3dAttribute BRICK_STATSYM file
	df = 284

	# get the critical t statistic value
	crit_t = scipy.stats.t.isf( conf[ "ana" ][ "loc_p" ], df )

	for hemi in [ "lh", "rh" ]:

		loc_stat_file = "%s_%s_reml-full.niml.dset[%s]" % ( paths[ "loc" ][ "glm" ],
		                                                    hemi,
		                                                    loc_stat_bricks
		                                                  )

		roi_file = "%s_%s-full.niml.dset" % ( paths[ "rois" ][ "dset" ], hemi )

		for ( roi_name, roi_val ) in loc_rois:

			roi_stat_file = "%s_%s_%s.txt" % ( paths[ "loc" ][ "roi_stat" ],
			                                   roi_name,
			                                   hemi
			                                 )

			xtr_cmd = [ "3dmaskdump",
			            "-mask", roi_file,
			            "-mrange", roi_val, roi_val,
			            "-o", roi_stat_file,
			            "-index",  # save the node index
			            "-noijk",  # j and k are just zeros
			            loc_stat_file
			          ]

			# 3dmaskdump won't overwrite, so need to manually remove any prior file
			if os.path.exists( roi_stat_file ):
				os.remove( roi_stat_file )

			fmri_tools.utils.run_cmd( xtr_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

			roi_stat = np.loadtxt( roi_stat_file )


			# add a column of zeros that will hold the parcel number
			roi_stat = np.hstack( ( roi_stat,
			                        np.zeros( ( roi_stat.shape[ 0 ], 1 ) )
			                      )
			                    )

			# index 0 is the node index
			lvf = [ roi_stat[ :, 1 ] > crit_t,
			        np.abs( roi_stat[ :, 1 ] ) <= crit_t,
			        roi_stat[ :, 1 ] < -crit_t
			      ]

			rvf = [ roi_stat[ :, 2 ] > crit_t,
			        np.abs( roi_stat[ :, 2 ] ) <= crit_t,
			        roi_stat[ :, 2 ] < -crit_t
			      ]

			i_parcel = 1

			for i_lvf in xrange( 3 ):
				for i_rvf in xrange( 3 ):

					in_parcel = np.logical_and( lvf[ i_lvf ], rvf[ i_rvf ] )

					assert( np.all( roi_stat[ in_parcel, 3 ] == 0 ) )

					roi_stat[ in_parcel, 3 ] = i_parcel

					i_parcel += 1

			assert( np.logical_not( np.any( roi_stat[ :, 3 ] == 0 ) ) )

			np.savetxt( roi_stat_file, roi_stat )

			roi_parc_file = "%s_%s_%s-full.niml.dset" % ( paths[ "loc" ][ "roi_parc" ],
			                                              roi_name,
			                                              hemi
			                                            )

			# now to convert it to a full niml dataset
			conv_cmd = [ "ConvertDset",
			             "-i_1D",
			             "-input", "%s[3]" % roi_stat_file,
			             "-o_niml",
			             "-prefix", roi_parc_file,
			             "-overwrite",
			             "-node_index_1D", "%s[0]" % roi_stat_file,
			             "-pad_to_node", "%d" % conf[ "subj" ][ "node_k" ][ hemi ]
			           ]

			fmri_tools.utils.run_cmd( conv_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

	for ( roi_name, _ ) in loc_rois:

		parcel_file = "%s_%s.txt" % ( paths[ "loc" ][ "parc" ], roi_name )

		parcel_count = np.zeros( ( 3, 3 ) )

		for hemi in [ "lh", "rh" ]:

			roi_stat_file = "%s_%s_%s.txt" % ( paths[ "loc" ][ "roi_stat" ],
			                                   roi_name,
			                                   hemi
			                                 )

			roi_stat = np.loadtxt( roi_stat_file )

			i_parcel = 1

			for i_lvf in xrange( 3 ):
				for i_rvf in xrange( 3 ):

					parcel_count[ i_lvf, i_rvf ] += np.sum( roi_stat[ :, 3 ] == i_parcel )

					i_parcel += 1

		np.savetxt( parcel_file, parcel_count, "%d" )


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
		loc_mask_file = "%s_%s-full.niml.dset" % ( paths[ "loc" ][ "mask" ],
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

			try:
				fmri_tools.utils.run_cmd( xtr_cmd,
				                          env = fmri_tools.utils.get_env(),
				                          log_path = paths[ "summ" ][ "log_file" ]
				                        )
			except:
				with file( roi_psc_file, "a" ):
					os.utime( roi_psc_file, None )

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
