"""
Set of routines to analyse single-subject fMRI data for the natural scenes
aperture fMRI experiment.
"""

from __future__ import division

import numpy as np


def extract_blocks( paths, conf ):
	"""Parse the average VTC for each ROI into the individual blocks.

	Parameters
	-----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'ns_aperture.config'.
	conf : dict
		Experiment configuration, as returned by 'get_conf' in
		'ns_aperture.config'.

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


def boot_cond_diff( paths, conf ):
	"""Compute the average and bootstrapped difference between the experiment
	   conditions.

	Parameters
	-----------
	paths : dict of strings
		Subject path structure, as returned by 'get_subj_paths' in
		'ns_aperture.config'.
	conf : dict
		Experiment configuration, as returned by 'get_conf' in
		'ns_aperture.config'.

	"""

	n_boot = 1000

	for roi_name in conf[ "ana" ][ "rois" ]:

		blk = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "block" ], roi_name ) )

		cond_a = blk[ blk[ :, 1 ] == 0, 0 ]
		cond_b = blk[ blk[ :, 1 ] == 1, 0 ]

		assert( len( cond_a ) == len( cond_b ) )

		n_samp = len( cond_a )

		diff_boot = np.zeros( ( n_boot ) )

		diff = np.mean( cond_a ) - np.mean( cond_b )

		for i_boot in xrange( n_boot ):

			cond_a_boot = cond_a[ np.random.randint( 0, n_samp, n_samp ) ]
			cond_b_boot = cond_b[ np.random.randint( 0, n_samp, n_samp ) ]

			diff_boot[ i_boot ] = np.mean( cond_a_boot ) - np.mean( cond_b_boot )

		cond_diff = np.hstack( ( diff,
		                         np.percentile( diff_boot, ( 2.5, 97.5 ) )
		                       )
		                     )

		np.save( "%s-%s.npy" % ( paths[ "ana" ][ "cond_diff" ], roi_name ),
		         cond_diff
		       )
