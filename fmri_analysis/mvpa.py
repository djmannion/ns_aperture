"""
Set of routines to run MVPA on the natural scenes
apertures fMRI experiment.
"""

from __future__ import division

import os, os.path
import itertools

import numpy as np
import scipy.stats

import fmri_tools.utils

import ns_aperture.fmri.exp


def filt( paths, conf ):
	"""Filters the timecourses for MVPA analysis"""

	i_surf_files = np.array( conf[ "subj" ][ "exp_runs" ] ).astype( "int" ) - 1

	surf_files = [ paths[ "func" ][ "surf_files" ][ i_surf ]
	               for i_surf in i_surf_files
	             ]

	filt_files = [ paths[ "mvpa" ][ "filt_files" ][ i_surf ]
	               for i_surf in i_surf_files
	             ]

	start_dir = os.getcwd()

	filt_dir = os.path.split( filt_files[ 0 ] )[ 0 ]

	os.chdir( filt_dir )

	for hemi in [ "lh", "rh" ]:

		for ( surf_file, filt_file ) in zip( surf_files, filt_files ):

			filt_cmd = [ "3dDetrend",
			             "-prefix", "%s_%s.niml.dset" % ( filt_file, hemi ),
			             "-polort", "%d" % ( conf[ "ana" ][ "poly_ord" ] - 1 ),
			             "-overwrite",
			             "%s_%s.niml.dset" % ( surf_file, hemi )
			           ]

			fmri_tools.utils.run_cmd( filt_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

			# convert to full
			full_filt_file = "%s_%s-full.niml.dset" % ( filt_file, hemi )
			pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

			fmri_tools.utils.sparse_to_full( "%s_%s.niml.dset" % ( filt_file, hemi ),
			                                 full_filt_file,
			                                 pad_node = pad_node,
			                                 log_path = paths[ "summ" ][ "log_file" ],
			                                 overwrite = True
			                               )

	os.chdir( start_dir )


def data_xtr( paths, conf ):
	"""Extract the timecourses for each node of each ROI"""

	os.chdir( paths[ "mvpa" ][ "base_dir" ] )

	# we use this in checks
	n_vols_per_run = int( conf[ "exp" ][ "run_len_s" ] /
	                      conf[ "acq" ][ "tr_s" ]
	                    )

	for ( roi_name, _ ) in conf[ "ana" ][ "rois" ]:

		data = []
		node_info = []

		for ( i_hemi, hemi ) in enumerate( [ "lh", "rh" ] ):

			# the mask file contains the parcel value for each node in the ROI
			mask_file = "%s_%s_%s-full.niml.dset" % ( paths[ "loc" ][ "roi_parc" ],
			                                          roi_name,
			                                          hemi
			                                        )

			run_data = []

			for ( i_run, filt_file ) in enumerate( paths[ "mvpa" ][ "filt_files" ] ):

				# input data is the (full) filtered timeseries for this run
				data_file = "%s_%s-full.niml.dset" % ( filt_file, hemi )

				# we only need this file temporarily
				out_file = os.path.join( paths[ "mvpa" ][ "base_dir" ], "temp.txt" )

				if os.path.exists( out_file ):
					os.remove( out_file )

				# use the ROI file to mask the input dataset
				xtr_cmd = [ "3dmaskdump",
				            "-mask", mask_file,
				            "-index",  # include the node index in output
				            "-noijk",
				            "-o", out_file,
				            data_file,
				            mask_file  # include the ROI parcel value in output
				          ]

				fmri_tools.utils.run_cmd( xtr_cmd,
				                          env = fmri_tools.utils.get_env(),
				                          log_path = paths[ "summ" ][ "log_file" ]
				                        )

				# the resulting data should be ( nodes, time + node index + parcel
				# value )
				curr_data = np.loadtxt( out_file )

				# check this claim
				assert( curr_data.shape[ 1 ] == ( n_vols_per_run + 2 ) )

				# clean up the unwanted output
				if os.path.exists( out_file ):
					os.remove( out_file )

				# append the timeseries data (ie. not including the node index or the
				# parcel value)
				run_data.append( curr_data[ :, 1:-1 ] )

				# we also want to know the node info, but since this is the same for
				# all runs we need only compute it once
				if i_run == 0:
					# we'd like to know which hemisphere it came from
					hemi_info = np.ones( ( curr_data.shape[ 0 ], 1 ) ) * i_hemi
					# and the node index and parcel value
					node_info.append( np.hstack( ( curr_data[ :, [ 0, -1 ] ], hemi_info ) ) )

			# run_data is a 10-item (runs) list of ( nodes x time ) arrays
			# converting it to array makes it ( run x nodes x time ) array
			run_data = np.array( run_data )

			# check shape is as expected
			assert( run_data.shape[ 0 ] == conf[ "subj" ][ "n_exp_runs" ] )
			assert( run_data.shape[ 2 ] == n_vols_per_run )

			data.append( np.array( run_data ) )

		# data is a 2-item (hemisphere) list of ( runs x nodes x time) arrays
		# want to join data along the nodes dimension; pass this axis (1) to
		# concatenate (which doesnt count the list dimension)
		# data becomes runs x nodes x time
		data = np.concatenate( data, axis = 1 )

		assert( data.shape[ 0 ] == conf[ "subj" ][ "n_exp_runs" ] )
		assert( data.shape[ 2 ] == n_vols_per_run )

		# node_info is a 2-item (hemisphere) list of nodes x 3 arrays
		# node_info becomes nodes (lh then rh) x 3 ( node index, parcel index,
		# hemi )
		node_info = np.vstack( node_info )

		# a region type (stimulated or 'blank') is defined by its parcel(s)
		for ( i_parc, parc_ind ) in enumerate( conf[ "ana" ][ "i_parc" ] ):

			# find the indices of nodes that correspond to any of the parcel indices
			# in this region type
			i_parc_nodes = np.any( [ node_info[ :, 1 ] == k
			                         for k in parc_ind
			                       ], axis = 0
			                     )

			# pull out the data and info from these nodes
			parc_data = data[ :, i_parc_nodes, : ]
			parc_info = node_info[ i_parc_nodes, : ]

			# ... and save
			parc_data_file = "%s_%s_%s.npy" % ( paths[ "mvpa" ][ "parc_data" ],
			                                    roi_name,
			                                    conf[ "ana" ][ "parc_lbl" ][ i_parc ]
			                                  )

			np.save( parc_data_file, parc_data )

			parc_info_file = "%s_%s_%s.npy" % ( paths[ "mvpa" ][ "parc_info" ],
			                                    roi_name,
			                                    conf[ "ana" ][ "parc_lbl" ][ i_parc ]
			                                  )

			np.save( parc_info_file, parc_info )



def info( paths, conf ):
	"""Extract the sequence information relevant to MVPA analysis"""

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

				# not the first or last block, which we want to discard
				if curr_block_num != 1 and curr_block_num != conf[ "exp" ][ "n_blocks" ]:

					cond = run_seq[ i_evt, seq_info[ "block_type" ] ]

					if cond == 0:
						block_type = 1
					else:
						block_type = -1

					# -1 to change from one-base to zero-base, -1 to adjust for
					# discarding the first block
					i_block = curr_block_num - 1 - 1

					run_info[ i_run, i_block, 1 ] = block_type

					time_s = run_seq[ i_evt, seq_info[ "time_s" ] ]

					time_vol = time_s / conf[ "acq" ][ "tr_s" ]

					# correct for the HRF delay
					time_adj = time_vol + conf[ "ana" ][ "hrf_corr_vol" ]

					run_info[ i_run, i_block, 0 ] = time_adj

	assert( np.sum( np.isnan( run_info ) ) == 0 )

	np.save( paths[ "mvpa" ][ "blk_data_info" ], run_info )


def blk_data_xtr( paths, conf ):
	"""Extracts blocked and z-scored data"""

	n_vol_per_block = int( conf[ "exp" ][ "block_len_s" ] /
	                       conf[ "acq" ][ "tr_s" ]
	                     )

	# blk_data_info is runs x blocks x details ( 0 = onset, 1 = cond )
	blk_data_info = np.load( paths[ "mvpa" ][ "blk_data_info" ] )

	for ( ( roi_name, _ ), parc_lbl ) in itertools.product( conf[ "ana" ][ "rois" ],
	                                                        conf[ "ana" ][ "parc_lbl" ]
	                                                      ):

		parc_data_file = "%s_%s_%s.npy" % ( paths[ "mvpa" ][ "parc_data" ],
		                                    roi_name,
		                                    parc_lbl
		                                  )

		# parc_data is runs x nodes x time
		parc_data = np.load( parc_data_file )

		( n_runs, n_nodes, _ ) = parc_data.shape
		( _, n_blocks, _ ) = blk_data_info.shape

		blk_data = np.empty( ( n_runs, n_blocks, n_nodes ) )
		blk_data.fill( np.NAN )

		for i_run in xrange( n_runs ):
			for i_block in xrange( n_blocks ):

				# onset volume for this block, HRF corrected
				i_block_start = blk_data_info[ i_run, i_block, 0 ]

				# get the volume range for this block
				i_vol_range = np.arange( i_block_start,
				                         i_block_start + n_vol_per_block,
				                       ).astype( "int" )

				# blk_tc becomes nodes x time
				blk_tc = parc_data[ i_run, ... ][ :, i_vol_range ]

				# just verify that is indeed the case
				assert( blk_tc.shape == ( n_nodes, len( i_vol_range ) ) )

				# ... because we want to average over the vol dimension
				blk_mean = np.mean( blk_tc, axis = 1 )

				# average over block timepoints
				blk_data[ i_run, i_block, : ] = blk_mean

			# z-score all the blocks in this run, for each node
			blk_data[ i_run, :, : ] = scipy.stats.zscore( blk_data[ i_run, :, : ],
			                                              axis = 0
			                                            )

		assert( np.sum( np.isnan( blk_data ) ) == 0 )

		blk_file = "%s_%s_%s.npy" % ( paths[ "mvpa" ][ "blk_data" ],
		                              roi_name,
		                              parc_lbl
		                            )

		np.save( blk_file, blk_data )


def rfe_classify( paths, conf ):
	"""Runs recursive feature elimination classification"""

	# load the info for each block, which doesn't depend on hemisphere
	# blk_data_info is runs x blocks x details ( 0 = onset, 1 = cond )
	blk_data_info = np.load( paths[ "mvpa" ][ "blk_data_info" ] )

	log_path = paths[ "mvpa" ][ "log_file" ]

	# run *indices* providing the training data for each fold
	train_info = [ [ 2, 3, 4, 5, 6, 7, 8, 9 ],
	               [ 0, 1, 4, 5, 6, 7, 8, 9 ],
	               [ 0, 1, 2, 3, 6, 7, 8, 9 ],
	               [ 0, 1, 2, 3, 4, 5, 8, 9 ],
	               [ 0, 1, 2, 3, 4, 5, 6, 7 ]
	             ]

	# ... and test data
	test_info = [ [ 0, 1 ],
	              [ 2, 3 ],
	              [ 4, 5 ],
	              [ 6, 7 ],
	              [ 8, 9 ]
	            ]

	for ( ( roi_name, _ ), parc_lbl ) in itertools.product( conf[ "ana" ][ "rois" ],
	                                                        conf[ "ana" ][ "parc_lbl" ]
	                                                      ):

		roi_parc_dir = os.path.join( paths[ "mvpa" ][ "rfe_base_dir" ],
		                             "%s_%s" % ( roi_name, parc_lbl )
		                           )

		blk_data_file = "%s_%s_%s.npy" % ( paths[ "mvpa" ][ "blk_data" ],
		                                   roi_name,
		                                   parc_lbl
		                                 )

		# blk_data is runs x blocks x nodes
		blk_data = np.load( blk_data_file )

		( _, n_blocks, n_nodes ) = blk_data.shape

		parc_info_file = "%s_%s_%s.npy" % ( paths[ "mvpa" ][ "parc_info" ],
		                                    roi_name,
		                                    parc_lbl
		                                  )

		# parc_info is nodes x 3 ( node index, parcel value, hemi )
		parc_info = np.load( parc_info_file )

		# this holds the indices for all the nodes
		i_nodes = parc_info[ :, 0 ]
		i_hemis = parc_info[ :, 2 ]

		assert( parc_info.shape[ 0 ] == n_nodes )

		svm_runs = zip( train_info, test_info )

		n_folds = len( svm_runs )

		# this will store our accuracy for each fold and RFE level
		acc = np.empty( ( n_folds, conf[ "ana" ][ "rfe_levels" ] ) )
		acc.fill( np.NAN )

		for ( i_fold, ( train_runs, test_runs ) ) in enumerate( svm_runs ):

			fold_dir = os.path.join( roi_parc_dir, "fold%02d" % ( i_fold + 1 ) )

			# check that there is no contamination of train and test sets
			assert( set.isdisjoint( set( train_runs ), set( test_runs ) ) )

			train_data = blk_data[ train_runs, :, : ]

			assert( train_data.shape == ( len( train_runs ), n_blocks, n_nodes ) )

			train_data_info = blk_data_info[ train_runs, :, 1 ]

			assert( train_data_info.shape == ( len( train_runs ), n_blocks ) )

			( n_train_runs, _, n_train_nodes ) = train_data.shape

			# make a copy of the node indices to use for this fold
			i_fold_nodes = i_nodes
			i_fold_hemis = i_hemis

			# how many splits to do per RFE level
			n_splits = int( n_train_runs / 2 )

			test_data = blk_data[ test_runs, :, : ]
			test_data_info = blk_data_info[ test_runs, :, 1 ]

			for i_level in xrange( conf[ "ana" ][ "rfe_levels" ] ):

				level_dir = os.path.join( fold_dir, "rfe%02d" % ( i_level + 1 ) )

				split_wgt = np.empty( ( n_train_nodes, n_splits ) )
				split_wgt.fill( np.NAN )

				for i_split in xrange( n_splits ):

					split_dir = os.path.join( level_dir, "split%02d" % ( i_split + 1 ) )

					# extract training data for this split
					i_split_runs = np.delete( np.arange( n_train_runs ),
					                          [ i_split * 2, i_split * 2 + 1 ]
					                        )

					split_data = train_data[ i_split_runs, :, : ]
					split_data_info = train_data_info[ i_split_runs, : ]

					split_wgt[ :, i_split ] = _get_split_weights( split_dir,
					                                              split_data,
					                                              split_data_info,
					                                              log_path
					                                            )

				# average weights over splits
				wgt = np.mean( split_wgt, axis = 1 )

				# get the indices that would sort |weights| (in ascending order )
				i_wgt = np.argsort( np.abs( wgt ) )

				# cull lowest performing weights
				n_to_cull = np.round( n_train_nodes * conf[ "ana" ][ "rfe_cull_p" ] )

				# update
				train_data = train_data[ :, :, i_wgt[ n_to_cull: ] ]
				n_train_nodes = train_data.shape[ -1 ]

				test_data = test_data[ :, :, i_wgt[ n_to_cull: ] ]
				n_test_nodes = test_data.shape[ -1 ]

				assert( n_train_nodes == n_test_nodes )

				# cull from the list of node indices
				i_fold_nodes = i_fold_nodes[ i_wgt[ n_to_cull: ] ]
				assert( len( i_fold_nodes ) == n_train_nodes )

				i_fold_hemis = i_fold_hemis[ i_wgt[ n_to_cull: ] ]
				assert( len( i_fold_hemis ) == n_train_nodes )

				# run the actual train / test
				[ lvl_acc, lvl_wgt ] = _run_level( level_dir,
				                                   train_data,
				                                   train_data_info,
				                                   test_data,
				                                   test_data_info,
				                                   log_path
				                                 )

				acc[ i_fold, i_level ] = lvl_acc

				node_and_wgt_info = np.column_stack( ( i_fold_nodes,
				                                       i_fold_hemis,
				                                       lvl_wgt
				                                     )
				                                   )

				node_and_wgt_file = os.path.join( level_dir,
				                                  "%s_%s_%s_fold_%02d.txt" % (
				                                    paths[ "mvpa" ]["node_wgt_base" ],
				                                    roi_name,
				                                    parc_lbl,
				                                    ( i_fold + 1 )
				                                  )
				                                )

				np.savetxt( node_and_wgt_file, node_and_wgt_info )

		acc_file = os.path.join( roi_parc_dir,
		                         "%s_%s_%s.txt" % ( paths[ "mvpa" ][ "acc_base" ],
		                                            roi_name,
		                                            parc_lbl
		                                          )
		                       )

		np.savetxt( acc_file, acc )


def _run_level( level_dir,
                train_data,
                train_data_info,
                test_data,
                test_data_info,
                log_path
              ):
	"""Run a train-and-test procedure for an RFE level"""

	# learn from all train data
	train_str = _get_svm_str( train_data, train_data_info )

	train_path = os.path.join( level_dir, "train.txt" )

	with open( train_path, "w" ) as train_file:
		train_file.write( train_str )

	model_path = os.path.join( level_dir, "model.txt" )
	weight_path = os.path.join( level_dir, "weight.txt" )

	train_cmd = [ "svm_learn",
	              train_path,
	              model_path
	            ]

	fmri_tools.utils.run_cmd( train_cmd,
	                          env = fmri_tools.utils.get_env(),
	                          log_path = log_path
	                        )

	# extract weights
	wgt_cmd = [ "svm2weight.pl",
	            model_path
	          ]

	wgt = fmri_tools.utils.run_cmd( wgt_cmd,
	                                env = fmri_tools.utils.get_env(),
	                                log_path = "return"
	                              )

	wgt = np.array( wgt.splitlines() ).astype( "float" )

	np.savetxt( weight_path, wgt )

	# apply to test data
	test_str = _get_svm_str( test_data, test_data_info )

	test_path = os.path.join( level_dir, "test.txt" )

	with open( test_path, "w" ) as test_file:
		test_file.write( test_str )

	pred_path = os.path.join( level_dir, "pred.txt" )

	test_cmd = [ "svm_classify",
	              test_path,
	              model_path,
	              pred_path
	            ]

	fmri_tools.utils.run_cmd( test_cmd,
	                          env = fmri_tools.utils.get_env(),
	                          log_path = log_path
	                        )

	pred = np.loadtxt( pred_path )

	true_conds = []

	for i_run in xrange( test_data_info.shape[ 0 ] ):
		for i_block in xrange( test_data_info.shape[ 1 ] ):

			ex_cond = test_data_info[ i_run, i_block ]

			true_conds.append( ex_cond )

	acc = np.sum( true_conds == np.sign( pred ) ) / len( pred ) * 100

	return [ acc, wgt ]


def _get_split_weights( split_dir, split_data, split_data_info, log_path ):
	"""Returns the weights from a given RFE split"""

	# learn SVM
	train_str = _get_svm_str( split_data, split_data_info )

	train_path = os.path.join( split_dir, "train.txt" )

	with open( train_path, "w" ) as train_file:
		train_file.write( train_str )

	model_path = os.path.join( split_dir, "model.txt" )
	weight_path = os.path.join( split_dir, "weight.txt" )

	train_cmd = [ "svm_learn",
	              train_path,
	              model_path
	            ]

	fmri_tools.utils.run_cmd( train_cmd,
	                          env = fmri_tools.utils.get_env(),
	                          log_path = log_path
	                        )

	# extract weights
	wgt_cmd = [ "svm2weight.pl",
	            model_path
	          ]

	wgt = fmri_tools.utils.run_cmd( wgt_cmd,
	                                env = fmri_tools.utils.get_env(),
	                                log_path = "return"
	                              )

	wgt = np.array( wgt.splitlines() ).astype( "float" )

	np.savetxt( weight_path, wgt )

	return wgt


def weight_maps( paths, conf ):
	"""Writes the weight maps to datasets"""

	for ( ( roi_name, _ ), parc_lbl ) in itertools.product( conf[ "ana" ][ "rois" ],
	                                                        conf[ "ana" ][ "parc_lbl" ]
	                                                      ):

		roi_parc_dir = os.path.join( paths[ "mvpa" ][ "rfe_base_dir" ],
		                             "%s_%s" % ( roi_name, parc_lbl )
		                           )

		os.chdir( roi_parc_dir )

		n_folds = 5

		for ( i_hemi, hemi ) in enumerate( [ "lh", "rh" ] ):

			hemi_dsets = []

			for i_fold in xrange( n_folds ):

				fold_dir = os.path.join( roi_parc_dir, "fold%02d" % ( i_fold + 1 ) )

				for i_level in xrange( conf[ "ana" ][ "rfe_levels" ] ):

					level_dir = os.path.join( fold_dir, "rfe%02d" % ( i_level + 1 ) )

					node_and_wgt_file = os.path.join( level_dir,
					                                  "%s_%s_%s_fold_%02d.txt" % (
					                                    paths[ "mvpa" ][ "node_wgt_base" ],
					                                    roi_name,
					                                    parc_lbl,
					                                    ( i_fold + 1 )
					                                  )
					                                )

					wgt_etc = np.loadtxt( node_and_wgt_file )

					wgt_etc_hemi = wgt_etc[ wgt_etc[ :, 1 ] == i_hemi, : ][ :, [ 0, 2 ] ]

					node_file = os.path.join( level_dir, "nodes_%s.1D" % hemi )
					np.savetxt( node_file, wgt_etc_hemi[ :, 0 ] )

					weight_file = os.path.join( level_dir, "weights_%s.1D" % hemi )
					np.savetxt( weight_file, wgt_etc_hemi[ :, 1 ] )

					dset_file = os.path.join( level_dir, "dset_%s-full.niml.dset" % hemi )

					conv_cmd = [ "ConvertDset",
					             "-i_1D",
					             "-input", weight_file,
					             "-node_index_1D", node_file,
					             "-o_niml",
					             "-prefix", dset_file,
					             "-pad_to_node", "%d" % conf[ "subj" ][ "node_k" ][ hemi ],
					             "-overwrite"
					           ]

					fmri_tools.utils.run_cmd( conv_cmd,
					                          env = fmri_tools.utils.get_env(),
					                          log_path = paths[ "summ" ][ "log_file" ]
					                        )

					hemi_dsets.append( dset_file )

			cat_file = os.path.join( roi_parc_dir,
			                         "cat_weights_%s-full.niml.dset" % hemi
			                       )

			# cat them all together
			cat_cmd = [ "3dTcat",
			            "-overwrite",
			            "-prefix", cat_file
			          ]

			cat_cmd.extend( hemi_dsets )

			fmri_tools.utils.run_cmd( cat_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

			mean_file = os.path.join( roi_parc_dir,
			                          "%s_%s_%s_%s-full.niml.dset" % (
			                            paths[ "mvpa" ][ "wgt_base" ],
			                            roi_name,
			                            parc_lbl,
			                            hemi
			                          )
			                        )

			mean_cmd = [ "3dMean",
			             "-overwrite",
			             "-prefix", mean_file
			           ]

			mean_cmd.extend( hemi_dsets )

			fmri_tools.utils.run_cmd( mean_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )


def _get_svm_str( blk_data, blk_data_info ):
	"""Returns a string in SVMlight format"""

	( n_runs, n_blocks, n_nodes ) = blk_data.shape

	assert( blk_data_info.shape[ 0 ] == n_runs )
	assert( blk_data_info.shape[ 1 ] == n_blocks )

	svm_str = ""

	for i_run in xrange( n_runs ):
		for i_block in xrange( n_blocks ):

			blk_cond = blk_data_info[ i_run, i_block ]

			svm_str += "%+d" % blk_cond

			for i_node in xrange( n_nodes ):

				svm_str += " %d:%.12f" % ( i_node + 1,
				                           blk_data[ i_run, i_block, i_node ]
				                         )

			svm_str += "\n"

	return svm_str


def weight_dist( paths, conf ):
	"""a"""

	# first, convert fovea ROI into a dataset
	for hemi in [ "lh", "rh" ]:

		dset_cmd = [ "ROI2dataset",
		             "-overwrite",
		             "-pad_to_node", "%d" % conf[ "subj" ][ "node_k" ][ hemi ],
		             "-prefix", "%s_%s-full.niml.dset" % ( paths[ "rois" ][ "fov_dset" ],
		                                                   hemi
		                                                 ),
		             "-input", os.path.join( paths[ "rois" ][ "base_dir" ],
		                                     "%s_fov.1D.roi" % hemi
		                                   )
		           ]

		fmri_tools.utils.run_cmd( dset_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		for ( roi_name, _ ) in conf[ "ana" ][ "rois" ]:

			os.chdir( os.path.split( paths[ "mvpa" ][ "wgt_fov" ] )[ 0 ] )

			roi_dir = os.path.join( paths[ "mvpa" ][ "rfe_base_dir" ],
		                          "%s_blnk" % roi_name
		                        )

			wgt_file = os.path.join( roi_dir,
			                         "%s_%s_blnk_%s-full.niml.dset" % (
			                           paths[ "mvpa" ][ "wgt_base" ],
			                           roi_name,
			                           hemi
			                         )
			                       )

			roi_file = "%s_%s-full.niml.dset" % ( paths[ "rois" ][ "fov_dset" ],
			                                      hemi
			                                    )

			base_cmd = [ "3dcalc",
			             "-overwrite",
			             "-a", wgt_file,
			             "-b", roi_file
			           ]

			fov_file = "%s_%s_%s-full.niml.dset" % ( paths[ "mvpa" ][ "wgt_fov" ],
			                                         roi_name,
			                                         hemi
			                                       )

			fov_cmd = list( base_cmd )
			fov_cmd.extend( [ "-expr", "a*ispositive(b)",
			                  "-prefix", fov_file
			                ]
			              )

			fmri_tools.utils.run_cmd( fov_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

			fov_txt = "%s_%s_%s.txt" % ( paths[ "mvpa" ][ "wgt_fov" ],
			                             roi_name,
			                             hemi
			                           )

			if os.path.exists( fov_txt ):
				os.remove( fov_txt )

			dump_cmd = [ "3dmaskdump",
			             "-mask", fov_file,
			             "-noijk",
			             "-o", fov_txt,
			             fov_file
			           ]

			fmri_tools.utils.run_cmd( dump_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

			ecc_file = "%s_%s_%s-full.niml.dset" % ( paths[ "mvpa" ][ "wgt_ecc" ],
			                                         roi_name,
			                                         hemi
			                                       )

			ecc_cmd = list( base_cmd )
			ecc_cmd.extend( [ "-expr", "a*not(ispositive(b))",
			                  "-prefix", ecc_file
			                ]
			              )

			fmri_tools.utils.run_cmd( ecc_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

			ecc_txt = "%s_%s_%s.txt" % ( paths[ "mvpa" ][ "wgt_ecc" ],
			                             roi_name,
			                             hemi
			                           )

			if os.path.exists( ecc_txt ):
				os.remove( ecc_txt )

			dump_cmd = [ "3dmaskdump",
			             "-mask", ecc_file,
			             "-noijk",
			             "-o", ecc_txt,
			             ecc_file
			           ]

			fmri_tools.utils.run_cmd( dump_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )


	fov_data = []
	ecc_data = []

	for ( i_roi, ( roi_name, _ ) ) in enumerate( conf[ "ana" ][ "rois" ] ):

		for ( i_hemi, hemi ) in enumerate( [ "lh", "rh" ] ):

			f = np.loadtxt( "%s_%s_%s.txt" % ( paths[ "mvpa" ][ "wgt_fov" ],
			                                   roi_name,
			                                   hemi
			                                 )
			              )[ :, np.newaxis ]

			h_f = np.ones( f.shape ) * i_hemi
			r_f = np.ones( f.shape ) * i_roi

			e = np.loadtxt( "%s_%s_%s.txt" % ( paths[ "mvpa" ][ "wgt_ecc" ],
			                                   roi_name,
			                                   hemi
			                                 )
			              )[ :, np.newaxis ]

			h_e = np.ones( e.shape ) * i_hemi
			r_e = np.ones( e.shape ) * i_roi

			fov_data.append( np.hstack( ( f, r_f, h_f ) ) )
			ecc_data.append( np.hstack( ( e, r_e, h_e ) ) )

	fov_data = np.vstack( fov_data )
	ecc_data = np.vstack( ecc_data )

	np.savetxt( "%s.txt" % paths[ "mvpa" ][ "wgt_fov" ], fov_data )
	np.savetxt( "%s.txt" % paths[ "mvpa" ][ "wgt_ecc" ], ecc_data )




