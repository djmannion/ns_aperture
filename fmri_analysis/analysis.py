"""
Set of routines to analyse single-subject fMRI data for the natural scenes
aperture fMRI experiment.
"""

from __future__ import division

import multiprocessing

import numpy as np


def loc_extract_blocks( paths, conf, subj_conf ):
	"""Parse the average VTC for each ROI into the individual blocks."""

	# volumes x runs x voxels
	vtc = np.load( "%s-gray.npy" % paths[ "ana_loc" ][ "vtc_file" ] )

	# blocks x runs x [ start, cond ]
	loc_design = np.load( paths[ "log" ][ "loc_design" ] )

	blk_data = np.empty( ( conf[ "exp" ][ "loc_n_valid_blocks" ] / 3,
	                       subj_conf[ "n_loc_runs" ],
	                       3,
	                       vtc.shape[ 2 ]
	                     )
	                   )
	blk_data.fill( np.NAN )

	block_cond_count = np.zeros( ( subj_conf[ "n_loc_runs" ], 3 ) )

	for i_run in xrange( blk_data.shape[ 1 ] ):

		run_vtc = vtc[ :, i_run, : ]

		run_design = loc_design[ :, i_run, : ]

		for i_blk in xrange( loc_design.shape[ 0 ] ):

			blk_cond = run_design[ i_blk, 1 ]
			i_blk_cond = block_cond_count[ i_run, blk_cond ]

			i_start = run_design[ i_blk, 0 ]
			i_end = i_start + conf[ "exp" ][ "n_vols_per_blk" ]

			blk_range = np.arange( i_start, i_end ).astype( "int" )

			blk_range += conf[ "ana" ][ "hrf_corr_vols" ]

			blk_mean = np.mean( run_vtc[ blk_range, : ], axis = 0 )

			blk_data[ i_blk_cond, i_run, blk_cond, : ] = blk_mean

			block_cond_count[ i_run, blk_cond ] += 1


	psc_data = np.empty( ( blk_data.shape[ 0 ] * 2,
	                       2,
	                       blk_data.shape[ -1 ]
	                     )
	                   )
	psc_data.fill( np.NAN )

	for i_run in xrange( blk_data.shape[ 1 ] ):

		baseline = np.mean( blk_data[ :, i_run, 0, : ], axis = 0 )

		for i_blk in xrange( blk_data.shape[ 0 ] ):

			i_blk_k = i_run * blk_data.shape[ 0 ] + i_blk

			for i_cond in xrange( 2 ):

				signal = blk_data[ i_blk, i_run, i_cond + 1, : ]

				psc_data[ i_blk_k, i_cond, : ] = ( signal - baseline ) / baseline * 100


	np.save( paths[ "ana_loc" ][ "block_psc_file" ], psc_data )


def loc_boot_analysis( paths, conf, subj_conf ):
	""""""

	# datapoints x conditions x voxels
	data = np.load( paths[ "ana_loc" ][ "block_psc_file" ] )

	n_boot = 100

	boot_data = np.empty( ( 3, data.shape[ -1 ], n_boot + 1 ) )

	n_samp = data.shape[ 0 ]

	cond_a = np.mean( data[ :, 0, : ], axis = 0 )
	cond_b = np.mean( data[ :, 1, : ], axis = 0 )

	boot_data[ 0, :, 0 ] = cond_a - cond_b
	boot_data[ 1, :, 0 ] = cond_a
	boot_data[ 2, :, 0 ] = cond_b

	n_processors = 16

	for i_boot in xrange( n_boot ):

		# open up our multiprocessor pool
		p = multiprocessing.Pool( processes = n_processors )

		b_data = p.map( _loc_bootstrap,
		                ( data[ :, :, i_voxel ]
		                  for i_voxel in xrange( data.shape[ -1 ] )
		                )
		              )

		boot_data[ :, :, i_boot + 1 ] = np.array( b_data ).T

		p.close()

	return boot_data


def _loc_bootstrap( data ):

	np.random.seed()

	# data is samples x conditions
	i_rand = np.random.randint( 0, data.shape[ 0 ], ( data.shape[ 0 ], 2 ) )

	cond_a = np.mean( data[ i_rand[ :, 0 ], 0 ] )
	cond_b = np.mean( data[ i_rand[ :, 1 ], 1 ] )

	return np.array( ( cond_a - cond_b, cond_a, cond_b ) )


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
