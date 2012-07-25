"""
Set of routines to analyse single-subject fMRI data for the natural scenes
aperture fMRI experiment.
"""

from __future__ import division

import os

import numpy as np

import fmri_tools.utils


def loc_glm( paths, conf ):
	"""Localiser GLM"""

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "loc_dir" ] )

	( lvf_file, rvf_file ) = paths[ "ana" ][ "loc_time_files" ]

	for hemi in [ "lh", "rh" ]:

		fit_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "loc_fits" ], hemi )
		glm_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "loc_glm" ], hemi )
		beta_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "loc_beta" ], hemi )

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		glm_cmd.extend( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                  for surf_file in paths[ "func_loc" ][ "surf_files" ]
		                ]
		              )

		glm_cmd.extend( [ "-force_TR", "%.3f" % conf[ "acq" ][ "tr_s" ],
		                  "-polort", "4",
		                  "-num_stimts", "2",
		                  "-stim_label", "1", "lvf_ON",
		                  "-stim_times", "1", lvf_file, "SPMG1(16)",
		                  "-stim_label", "2", "rvf_ON",
		                  "-stim_times", "2", rvf_file, "SPMG1(16)",
		                  "-local_times",
		                  "-gltsym", "SYM: +lvf_ON -rvf_ON",
		                  "-glt_label", "1", "LVFgtRVF",
		                  "-gltsym", "SYM: +lvf_ON \ +rvf_ON",
		                  "-glt_label", "2", "ALLf",
		                  "-xjpeg", "loc_design.png",
		                  "-x1D", "loc_design",
		                  "-jobs", "16",
		                  "-fitts", fit_file,
		                  "-bucket", glm_file,
		                  "-cbucket", beta_file,
		                  "-tout",
		                  "-fout",
		                  "-overwrite"
		                ]
		              )

		fmri_tools.utils.run_cmd( glm_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		fmri_tools.utils.sparse_to_full( glm_file,
		                                 "%s_%s-full" % ( paths[ "ana" ][ "loc_glm" ], hemi ),
		                                 pad_node = pad_node,
		                                 log_path = paths[ "summ" ][ "log_file" ],
		                                 overwrite = True
		                               )

		# run FDR correction
		fdr_cmd = [ "3dFDR",
		            "-input", glm_file,
		            "-prefix", "%s_%s.niml.dset" % ( paths[ "ana" ][ "loc_q" ], hemi ),
		            "-qval",
		            "-float",
		            "-overwrite"
		          ]

		fmri_tools.utils.run_cmd( fdr_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		q_file = "%s_%s" % ( paths[ "ana" ][ "loc_q" ], hemi )


		# convert the FDR to full
		fmri_tools.utils.sparse_to_full( "%s.niml.dset" % q_file,
		                                 "%s-full" % q_file,
		                                 pad_node = pad_node,
		                                 log_path = paths[ "summ" ][ "log_file" ],
		                                 overwrite = True
		                               )

		# use it to mask the ROIs
		mask_cmd = [ "3dcalc",
		             "-a", "%s-full.niml.dset[6]" % q_file,
		             "-b", "%s_%s-full.niml.dset" % (
		                     paths[ "ana" ][ "roi_dset" ],
		                     hemi ),
		             "-expr", "within( a, 0, %.5f ) * b" % conf[ "ana" ][ "q_thr" ],
		             "-prefix", "%s_%s-full.niml.dset" % (
		                          paths[ "ana" ][ "roi_mask" ],
		                          hemi ),
		             "-overwrite"
		           ]

		fmri_tools.utils.run_cmd( mask_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )


	os.chdir( start_dir )


def exp_glm( paths, conf ):
	"""Experiment GLM"""

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "exp_dir" ] )

	exp_file = paths[ "ana" ][ "exp_time_file" ]

	for hemi in [ "lh", "rh" ]:

		fit_file = "%s_%s" % ( paths[ "ana" ][ "exp_fits" ], hemi )
		glm_file = "%s_%s" % ( paths[ "ana" ][ "exp_glm" ], hemi )
		beta_file = "%s_%s" % ( paths[ "ana" ][ "exp_beta" ], hemi )

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		glm_cmd.extend( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                  for surf_file in paths[ "func_exp" ][ "surf_files" ]
		                ]
		              )

		glm_cmd.extend( [ "-force_TR", "%.3f" % conf[ "acq" ][ "tr_s" ],
		                  "-polort", "4",
		                  "-num_stimts", "3",
		                  "-stim_label", "1", "coh",
		                  "-stim_times", "1", exp_file, "SPMG1(16)",
		                  "-stim_label", "2", "noint-Aon",
		                  "-stim_file", "2", paths[ "log" ][ "reg_A" ],
		                  "-stim_label", "3", "noint-Bon",
		                  "-stim_file", "3", paths[ "log" ][ "reg_B" ],
		                  "-local_times",
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-jobs", "16",
		                  "-fitts", "%s.niml.dset" % fit_file,
		                  "-bucket", "%s.niml.dset" % glm_file,
		                  "-cbucket", "%s.niml.dset" % beta_file,
		                  "-vout",
		                  "-tout",
		                  "-overwrite",
		                ]
		              )

		fmri_tools.utils.run_cmd( glm_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		# convert the output to full
		for out_file in [ fit_file, glm_file, beta_file ]:

			fmri_tools.utils.sparse_to_full( "%s.niml.dset" % out_file,
			                                 "%s-full" % out_file,
			                                 pad_node = pad_node,
			                                 log_path = paths[ "summ" ][ "log_file" ],
			                                 overwrite = True
			                               )


	os.chdir( start_dir )


def beta_to_psc( paths, conf ):
	"""a"""

	exp_type = [ "exp", "loc" ]

	i_betas = [ "50", "10,11" ]

	start_dir = os.getcwd()

	for ( i_et, et ) in enumerate( exp_type ):

		dir_key = "%s_dir" % et

		os.chdir( paths[ "ana" ][ dir_key ] )

		for hemi in [ "lh", "rh" ]:

			beta_key = "%s_beta" % et
			beta_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ beta_key ], hemi )

			mat_file = os.path.join( paths[ "ana" ][ dir_key ],
		                           "%s_design.xmat.1D" % et
		                         )

			bltc_key = "%s_bltc" % et
			bltc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ bltc_key ], hemi )

			# generate an average baseline timecourse
			bl_cmd = [ "3dSynthesize",
			           "-cbucket", beta_file,
			           "-matrix", mat_file,
			           "-select", "baseline",
			           "-prefix", bltc_file,
			           "-overwrite"
			         ]

			fmri_tools.utils.run_cmd( bl_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

			bl_key = "%s_bl" % et
			bl_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ bl_key ], hemi )

			# average baseline across time
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

			psc_key = "%s_psc" % et
			psc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ psc_key ], hemi )

			# compute psc
			# from http://afni.nimh.nih.gov/sscc/gangc/TempNorm.html
			psc_cmd = [ "3dcalc",
			            "-fscale",
			            "-a", bl_file,
			            "-b", "%s[%s]" % ( beta_file, i_betas[ i_et ] ),
			            "-expr", "100 * b/a * step( 1 - abs( b/a ) )",
			            "-prefix", psc_file,
			            "-overwrite"
			          ]

			fmri_tools.utils.run_cmd( psc_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

			# convert to full
			full_psc_file = "%s_%s-full" % ( paths[ "ana" ][ psc_key ], hemi )

			pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

			fmri_tools.utils.sparse_to_full( psc_file,
			                                 full_psc_file,
			                                 pad_node = pad_node,
			                                 log_path = paths[ "summ" ][ "log_file" ],
			                                 overwrite = True
			                               )


def roi_xtr( paths, conf ):
	"""a"""

	for hemi in [ "lh", "rh" ]:

		# the *full* ROI file
		roi_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "roi_dset" ], hemi )

		for ( roi_name, roi_val ) in conf[ "ana" ][ "rois" ]:

			for ana_type in [ "exp", "loc" ]:

				for data_type in [ "psc", "stat" ]:

					if data_type == "psc":
						data_key = "%s_psc" % ana_type
					else:
						data_key = "%s_glm" % ana_type

					roi_key = "%s_roi_%s" % ( ana_type, data_type )

					data_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ data_key ],
					                                       hemi
					                                     )

					out_file = "%s_%s_%s.txt" % ( paths[ "ana" ][ roi_key ],
					                              roi_name,
					                              hemi
					                            )

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


def get_adj_tc( paths, conf ):
	"""UNFINISHED"""

	exp_type = [ "exp", "loc" ]

	i_betas = [ "50", "10,11" ]

	start_dir = os.getcwd()

	for ( i_et, et ) in enumerate( exp_type ):

		dir_key = "%s_dir" % et

		os.chdir( paths[ "ana" ][ dir_key ] )

		for hemi in [ "lh", "rh" ]:

			bltc_key = "%s_bltc" % et
			bltc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ bltc_key ], hemi )

			# generate an average baseline timecourse
			bl_cmd = [ "3dSynthesize",
			           "-cbucket", beta_file,
			           "-matrix", mat_file,
			           "-select", "baseline",
			           "-prefix", bltc_file,
			           "-overwrite"
			         ]

			fmri_tools.utils.run_cmd( bl_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

	os.chdir( start_dir )
