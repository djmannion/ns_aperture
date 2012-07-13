"""
Set of routines to analyse single-subject fMRI data for the natural scenes
aperture fMRI experiment.
"""

from __future__ import division

import multiprocessing

import numpy as np

import nipy


def loc_boot_analysis( paths, conf, subj_conf ):
	""""""

	# load the gray coords
	gray_coords = np.load( paths[ "roi" ][ "gray_coord_file" ] )

	n_voxels = gray_coords.shape[ 1 ]

	boot_data = np.load( paths[ "ana_loc" ][ "block_boot_file" ] )

	p_val = conf[ "ana" ][ "loc_p_thresh" ] * 100 / 2.0

	p_thresh = [ p_val, 100 - p_val ]

	# don't include the first index, which is the original data
	ci = np.array( np.percentile( boot_data[ :, :, 1: ],
	                              p_thresh,
	                              axis = 2
	                            )
	             )

	# ( L <> R, L > 0, R > 0, [ L | R ] > 0 )
	sig = np.empty( ( 4, n_voxels ) )
	sig.fill( np.NAN )

	# L <> R
	#  first test if the difference interval includes zero, then apply a NOT
	#  operation to return those outside the zero interval
	sig[ 0, : ] = np.logical_not( np.logical_and( ci[ 0, 0, : ] < 0,
	                                              0 < ci[ 1, 0, : ]
	                                            )
	                            )

	# gt 0 comparisons are signed, so can just check the lower bound to see if it
	# is greater than zero
	sig[ 1:3, : ] = ( ci[ 0, 1:, : ] > 0 )

	# combined is whether either L or R is greater than zero
	sig[ 3, : ] = np.logical_or( sig[ 1, : ], sig[ 2, : ] )

	np.save( paths[ "ana_loc" ][ "sig" ], sig )


def loc_write_img( paths, conf, subj_conf ):
	""""""

	# load the base image
	base_img = nipy.load_image( "%s.nii" % paths[ "summ" ][ "mean_file" ] )

	# initalise the volumes
	l_gt_r = np.zeros( ( base_img.shape ) )
	l_gt_z = l_gt_r.copy()
	r_gt_z = l_gt_r.copy()

	for i_voxel in xrange( gray_coords.shape[ 1 ] ):

		if sig[ 0, i_voxel ] == 1:
			l_gt_r[ gray_coords[ 0, i_voxel ],
			        gray_coords[ 1, i_voxel ],
			        gray_coords[ 2, i_voxel ]
			      ] = boot_data[ 0, i_voxel, 0 ]

		if sig[ 1, i_voxel ] == 1:
			l_gt_z[ gray_coords[ 0, i_voxel ],
			        gray_coords[ 1, i_voxel ],
			        gray_coords[ 2, i_voxel ]
			      ] = boot_data[ 1, i_voxel, 0 ]

		if sig[ 2, i_voxel ] == 1:
			r_gt_z[ gray_coords[ 0, i_voxel ],
			        gray_coords[ 1, i_voxel ],
			        gray_coords[ 2, i_voxel ]
			      ] = boot_data[ 2, i_voxel, 0 ]

	comb = np.logical_or( l_gt_z > 0, r_gt_z > 0 ).astype( "float" )

	l_gt_r_img = nipy.core.api.Image( l_gt_r, base_img.coordmap )
	nipy.save_image( l_gt_r_img, paths[ "ana_loc" ][ "l_gt_r_img" ] )

	l_gt_z_img = nipy.core.api.Image( l_gt_z, base_img.coordmap )
	nipy.save_image( l_gt_z_img, paths[ "ana_loc" ][ "l_gt_z_img" ] )

	r_gt_z_img = nipy.core.api.Image( r_gt_z, base_img.coordmap )
	nipy.save_image( r_gt_z_img, paths[ "ana_loc" ][ "r_gt_z_img" ] )

	comb_img = nipy.core.api.Image( comb, base_img.coordmap )
	nipy.save_image( comb_img, paths[ "ana_loc" ][ "comb_img" ] )



def loc_extract_blocks( paths, conf, subj_conf ):
	"""Parse the VTC for each voxel into the individual blocks."""

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


def loc_boot( paths, conf, subj_conf ):
	"""Run bootstrapping of the localiser data"""

	# datapoints x conditions x voxels
	data = np.load( paths[ "ana_loc" ][ "block_psc_file" ] )

	( n_samp, _, n_voxels ) = data.shape

	n_boot = conf[ "ana" ][ "n_boot" ]

	boot_data = np.empty( ( 3, n_voxels, n_boot + 1 ) )

	cond_a = np.mean( data[ :, 0, : ], axis = 0 )
	cond_b = np.mean( data[ :, 1, : ], axis = 0 )

	# initial, unbootstrapped data goes in slot 0
	boot_data[ 0, :, 0 ] = cond_a - cond_b
	boot_data[ 1, :, 0 ] = cond_a
	boot_data[ 2, :, 0 ] = cond_b

	n_processors = 16

	np.random.seed( subj_conf[ "boot_seed" ][ 1 ] )

	# generate a set of random indices for each condition and voxel
	i_rand = np.random.randint( 0,
	                            n_samp,
	                            ( n_boot, n_samp, 2, n_voxels )
	                          )

	for i_boot in xrange( n_boot ):

		# open up our multiprocessor pool
		p = multiprocessing.Pool( processes = n_processors )

		# distribute the voxels over processors
		b_data = p.map( _loc_bootstrap,
		                ( ( data[ :, :, i_voxel ],
		                    i_rand[ i_boot, :, :, i_voxel ]
		                  )
		                  for i_voxel in xrange( n_voxels )
		                )
		              )

		boot_data[ :, :, i_boot + 1 ] = np.array( b_data ).T

		p.close()

	np.save( paths[ "ana_loc" ][ "block_boot_file" ], boot_data )

	np.random.seed()


def _loc_bootstrap( args ):
	"""Helper function for localiser bootstrapping analysis."""

	( data, i_rand ) = args

	cond_a = np.mean( data[ i_rand[ :, 0 ], 0 ] )
	cond_b = np.mean( data[ i_rand[ :, 1 ], 1 ] )

	return np.array( ( cond_a - cond_b, cond_a, cond_b ) )


def exp_extract_blocks( paths, conf, subj_conf ):
	"""Parse the VTC for each voxel into the individual blocks."""

	# volumes x runs x voxels
	vtc = np.load( "%s-gray.npy" % paths[ "ana_exp" ][ "vtc_file" ] )

	# blocks x runs x [ start, cond ]
	exp_design = np.load( paths[ "log" ][ "design" ] )

	blk_data = np.empty( ( conf[ "exp" ][ "n_valid_blocks" ] / 2,
	                       subj_conf[ "n_runs" ],
	                       2,  # conditions
	                       vtc.shape[ 2 ]  # voxels
	                     )
	                   )
	blk_data.fill( np.NAN )

	block_cond_count = np.zeros( ( subj_conf[ "n_runs" ], 2 ) )

	for i_run in xrange( blk_data.shape[ 1 ] ):

		run_vtc = vtc[ :, i_run, : ]

		# convert to psc
		run_baseline = np.mean( run_vtc, axis = 0 )

		run_vtc = ( run_vtc - run_baseline ) / run_baseline * 100.0

		run_design = exp_design[ :, i_run, : ]

		for i_blk in xrange( exp_design.shape[ 0 ] ):

			blk_cond = run_design[ i_blk, 1 ]
			i_blk_cond = block_cond_count[ i_run, blk_cond ]

			i_start = run_design[ i_blk, 0 ]
			i_end = i_start + conf[ "exp" ][ "n_vols_per_blk" ]

			blk_range = np.arange( i_start, i_end ).astype( "int" )

			blk_range += conf[ "ana" ][ "hrf_corr_vols" ]

			blk_mean = np.mean( run_vtc[ blk_range, : ], axis = 0 )

			blk_data[ i_blk_cond, i_run, blk_cond, : ] = blk_mean

			block_cond_count[ i_run, blk_cond ] += 1

	# collapse over runs (column-major order)
	blk_data = np.reshape( blk_data,
	                       ( blk_data.shape[ 0 ] * blk_data.shape[ 1 ],
	                         blk_data.shape[ 2 ], blk_data.shape[ 3 ]
	                       ),
	                       order = "F"
	                     )

	np.save( paths[ "ana_exp" ][ "block_psc_file" ], blk_data )


def exp_boot( paths, conf, subj_conf ):
	"""Run bootstrapping of the experiment data"""

	# datapoints x conditions x voxels
	data = np.load( paths[ "ana_exp" ][ "block_psc_file" ] )

	( n_samp, _, n_voxels ) = data.shape

	n_boot = conf[ "ana" ][ "n_boot" ]

	boot_data = np.empty( ( n_voxels, n_boot + 1 ) )

	cond_a = np.mean( data[ :, 0, : ], axis = 0 )
	cond_b = np.mean( data[ :, 1, : ], axis = 0 )

	# initial, unbootstrapped data goes in slot 0
	boot_data[ :, 0 ] = cond_a - cond_b

	n_processors = 16

	np.random.seed( subj_conf[ "boot_seed" ][ 0 ] )

	# generate a set of random indices for each condition and voxel
	i_rand = np.random.randint( 0,
	                            n_samp,
	                            ( n_boot, n_samp, 2, n_voxels )
	                          )

	for i_boot in xrange( n_boot ):


		# open up our multiprocessor pool
		p = multiprocessing.Pool( processes = n_processors )

		# distribute the voxels over processors
		b_data = p.map( _exp_bootstrap,
		                ( ( data[ :, :, i_voxel ],
		                    i_rand[ i_boot, :, :, i_voxel ]
		                  )
		                  for i_voxel in xrange( n_voxels )
		                )
		              )

		boot_data[ :, i_boot + 1 ] = np.array( b_data )

		p.close()

	np.save( paths[ "ana_exp" ][ "block_boot_file" ], boot_data )

	np.random.seed()


def _exp_bootstrap( args ):
	"""Helper function for exp bootstrapping analysis."""

	( data, i_rand ) = args

	cond_a = np.mean( data[ i_rand[ :, 0 ], 0 ] )
	cond_b = np.mean( data[ i_rand[ :, 1 ], 1 ] )

	return ( cond_a - cond_b )


