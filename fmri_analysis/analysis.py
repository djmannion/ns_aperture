"""
Set of routines to analyse single-subject fMRI data for the natural scenes
aperture fMRI experiment.
"""

from __future__ import division

import os.path

import nipy
import numpy as np
import scipy.io

import fmri_tools.preproc, fmri_tools.utils


def extract_blocks( paths, conf, subj_conf ):
	"""Parse the average VTC for each ROI into the individual blocks.

	Parameters
	-----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'ns_aperture.config'.
	conf : dict
		Experiment configuration, as returned by 'get_conf' in
		'ns_aperture.config'.
	subj_conf : dict
		Subject configuration, as returned by 'get_subj_conf' in
		'ns_aperture.config', for this subject.

	"""

	design = np.load( paths[ "design" ][ "design" ] )

	for roi_name in conf[ "ana" ][ "rois" ]:

		vtc = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "vtc_avg" ],
		                               roi_name
		                             )
		             )

		data = []

		for i_run in xrange( vtc.shape[ 1 ] ):

			run_design = design[ :, i_run, : ]

			for i_blk in xrange( run_design.shape[ 0 ] ):

				blk_cond = run_design[ i_blk, 1 ]

				i_blk_start = run_design[ i_blk, 0 ]
				i_blk_end = i_blk_start + conf[ "exp" ][ "n_vols_per_blk" ]

				i_blk_range = np.arange( i_blk_start, i_blk_end ).astype( "int" )

				blk_mean = np.mean( vtc[ i_blk_range, i_run ] )

				data.append( [ blk_mean, blk_cond, i_blk, i_run ] )

		data = np.array( data )

		np.save( "%s-%s" % ( paths[ "ana" ][ "block" ], roi_name ),
		         data
		       )
