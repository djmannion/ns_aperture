"""
Set of routines to analyse single-subject fMRI data for the natural scenes
aperture fMRI experiment.
"""

from __future__ import division

import os

import numpy as np

import nipy

import fmri_tools.utils


def loc_glm( paths, conf ):
	"""Localiser GLM"""

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "loc_dir" ] )

	( lvf_file, rvf_file ) = paths[ "ana" ][ "loc_time_files" ]

	for hemi in [ "lh", "rh" ]:

		fit_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "loc_fits" ], hemi )
		glm_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "loc_glm" ], hemi )

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
		                  "-gltsym", "SYM: +lvf_ON",
		                  "-glt_label", "1", "LVF",
		                  "-gltsym", "SYM: +rvf_ON",
		                  "-glt_label", "2", "RVF",
		                  "-gltsym", "SYM: +lvf_ON -rvf_ON",
		                  "-glt_label", "3", "LVFgtRVF",
		                  "-xjpeg", "loc_design.png",
		                  "-x1D", "loc_design.txt",
		                  "-jobs", "16",
		                  "-fitts", fit_file,
		                  "-bucket", glm_file,
		                  "-overwrite"
		                ]
		              )

		fmri_tools.utils.run_cmd( glm_cmd,
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

		fit_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "exp_fits" ], hemi )
		glm_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "exp_glm" ], hemi )
		beta_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "exp_beta" ], hemi )

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		glm_cmd.extend( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                  for surf_file in paths[ "func_exp" ][ "surf_files" ]
		                ]
		              )

		glm_cmd.extend( [ "-force_TR", "%.3f" % conf[ "acq" ][ "tr_s" ],
		                  "-polort", "4",
		                  "-num_stimts", "1",
		                  "-stim_label", "1", "coh",
		                  "-stim_times", "1", exp_file, "SPMG1(16)",
		                  "-local_times",
		                  "-gltsym", "SYM: +coh",
		                  "-glt_label", "1", "coh",
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design.txt",
		                  "-jobs", "16",
		                  "-fitts", fit_file,
		                  "-bucket", glm_file,
		                  "-cbucket", beta_file,
		                  "-overwrite",
		                ]
		              )

		fmri_tools.utils.run_cmd( glm_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

	os.chdir( start_dir )
