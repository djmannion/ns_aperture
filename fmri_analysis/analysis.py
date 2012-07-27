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

		fit_file = "%s_%s" % ( paths[ "ana" ][ "loc_fits" ], hemi )
		glm_file = "%s_%s" % ( paths[ "ana" ][ "loc_glm" ], hemi )
		beta_file = "%s_%s" % ( paths[ "ana" ][ "loc_beta" ], hemi )

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
		                  "-glt_label", "2", "fALL",
		                  "-xjpeg", "loc_design.png",
		                  "-x1D", "loc_design",
		                  "-jobs", "16",
		                  "-fitts", fit_file,
		                  "-bucket", glm_file,
		                  "-cbucket", beta_file,
		                  "-tout",
		                  "-fout",
		                  "-overwrite",
		                  "-x1D_stop"
		                ]
		              )

		fmri_tools.utils.run_cmd( glm_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		reml_cmd = [ "3dREMLfit",
		             "-matrix", "loc_design.xmat.1D",
		             "-input", " ".join( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                                   for surf_file in paths[ "func_loc" ][ "surf_files" ]
		                                 ]
		                               ),
		             "-Rbeta", "%s_reml.niml.dset" % beta_file,
		             "-tout",
		             "-fout",
		             "-Rbuck", "%s_reml.niml.dset" % glm_file,
		             "-Rfitts", "%s_reml.niml.dset" % fit_file
		           ]

		fmri_tools.utils.run_cmd( reml_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		fmri_tools.utils.sparse_to_full( glm_file,
		                                 "%s_%s_reml-full" % ( paths[ "ana" ][ "loc_glm" ], hemi ),
		                                 pad_node = pad_node,
		                                 log_path = paths[ "summ" ][ "log_file" ],
		                                 overwrite = True
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
		                  "-fitts", "%s" % fit_file,
		                  "-bucket", "%s" % glm_file,
		                  "-cbucket", "%s" % beta_file,
		                  "-vout",
		                  "-tout",
		                  "-overwrite",
		                  "-x1D_stop"
		                ]
		              )

		fmri_tools.utils.run_cmd( glm_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		reml_cmd = [ "3dREMLfit",
		             "-matrix", "exp_design.xmat.1D",
		             "-input", " ".join( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                                   for surf_file in paths[ "func_exp" ][ "surf_files" ]
		                                 ]
		                               ),
		             "-Rbeta", "%s_reml.niml.dset" % beta_file,
		             "-tout",
		             "-Rbuck", "%s_reml.niml.dset" % glm_file,
		             "-Rfitts", "%s_reml.niml.dset" % fit_file
		           ]

		fmri_tools.utils.run_cmd( reml_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		# convert the output to full
		for out_file in [ fit_file, glm_file, beta_file ]:

			fmri_tools.utils.sparse_to_full( "%s_reml.niml.dset" % out_file,
			                                 "%s_reml-full" % out_file,
			                                 pad_node = pad_node,
			                                 log_path = paths[ "summ" ][ "log_file" ],
			                                 overwrite = True
			                               )

	os.chdir( start_dir )


def beta_to_psc( paths, conf ):
	"""Convert the GLM beta weights into units of percent signal change"""

	# want to convert both the experiment and the localiser data
	exp_type = [ "exp", "loc" ]

	# these are the indices into the beta files for the data we want to convert
	# they are checked below
	i_betas = [ "50", "10,11" ]

	start_dir = os.getcwd()

	# loop over the experiment / localiser experiment
	for ( i_et, et ) in enumerate( exp_type ):

		dir_key = "%s_dir" % et

		os.chdir( paths[ "ana" ][ dir_key ] )

		for hemi in [ "lh", "rh" ]:

			# dataset holding the beta weights
			beta_key = "%s_beta" % et
			beta_file = "%s_%s_reml.niml.dset" % ( paths[ "ana" ][ beta_key ], hemi )

			# design matrix file
			mat_file = os.path.join( paths[ "ana" ][ dir_key ],
			                         "%s_design.xmat.1D" % et
			                       )

			# baseline timecourse dataset, to write
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

			# baseline (point-estimate) dataset, to write
			bl_key = "%s_bl" % et
			bl_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ bl_key ], hemi )

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
			psc_key = "%s_psc" % et
			psc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ psc_key ], hemi )

			# the input beta file, with sub-brick selector
			beta_sel = "%s[%s]" % ( beta_file, i_betas[ i_et ] )

			# check that the label is as expected
			beta_label = fmri_tools.utils.get_dset_label( beta_sel )

			if et == "exp":
				assert( beta_label == "coh#0" )
			elif et == "loc":
				assert( beta_label == [ "lvf_ON#0", "rvf_ON#0" ] )

			# compute psc
			# from http://afni.nimh.nih.gov/sscc/gangc/TempNorm.html
			psc_cmd = [ "3dcalc",
			            "-fscale",
			            "-a", bl_file,
			            "-b", beta_sel,
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
	"""Extract PSC and statistics data from ROIs"""

	for hemi in [ "lh", "rh" ]:

		# the *full* ROI file
		roi_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "roi_dset" ], hemi )

		# iterate over all the ROIs
		for ( roi_name, roi_val ) in conf[ "ana" ][ "rois" ]:

			# ... and experiment types
			for ana_type in [ "exp", "loc" ]:

				# ... and what data we want to pull out
				for data_type in [ "psc", "stat" ]:

					if data_type == "psc":
						data_key = "%s_psc" % ana_type
						extra_desc = ""
					else:
						data_key = "%s_glm" % ana_type
						extra_desc = "_reml"

					roi_key = "%s_roi_%s" % ( ana_type, data_type )

					# our input dataset - either psc or glm
					data_file = "%s_%s%s-full.niml.dset" % ( paths[ "ana" ][ data_key ],
					                                         hemi,
					                                         extra_desc
					                                       )

					# our output dataset
					out_file = "%s_%s_%s.txt" % ( paths[ "ana" ][ roi_key ],
					                              roi_name,
					                              hemi
					                            )

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
