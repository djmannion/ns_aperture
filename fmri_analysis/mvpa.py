"""
Set of routines to run MVPA on the natural scenes
apertures fMRI experiment.
"""

from __future__ import division

import os, os.path
import tempfile

import numpy as np
import scipy.stats
import progressbar

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

	for hemi in [ "lh", "rh" ]:

		data = []

		# the *full* ROI file for this hemisphere
		roi_file = "%s_%s-full.niml.dset" % ( paths[ "rois" ][ "dset" ], hemi )

		for filt_file in paths[ "mvpa" ][ "filt_files" ]:

			# our input dataset
			data_file = "%s_%s-full.niml.dset" % ( filt_file, hemi )

			out_file = os.path.join( paths[ "mvpa" ][ "base_dir" ], "temp.txt" )

			if os.path.exists( out_file ):
				os.remove( out_file )

			# use the ROI file to mask the input dataset
			xtr_cmd = [ "3dmaskdump",
			            "-mask", roi_file,
			            "-index",  # include the node index in output
			            "-o", out_file,
			            data_file,
			            roi_file  # include the ROI value in output
			          ]

			fmri_tools.utils.run_cmd( xtr_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

			# append the output to this hemisphere's data
			data.append( np.loadtxt( out_file ) )

			# clean up the unwanted output
			if os.path.exists( out_file ):
				os.remove( out_file )

		# both the items in `data` are arrays with the same shape, so we can do a
		# straight conversion to a numpy array, which will be runs x nodes x time
		data = np.array( data )

		# pull put the node index (first entry) and roi index (last entry)
		data_info = data[ 0, :, [ 0, -1 ] ].astype( "int" )

		# ... and remove them the data array
		data = data[ :, :, 1:-1 ]

		# now we're ready to save
		data_file = "%s-%s.npy" % ( paths[ "mvpa" ][ "data" ], hemi )

		np.save( data_file, data )

		data_info_file = "%s-%s.npy" % ( paths[ "mvpa" ][ "data_info" ], hemi )

		np.save( data_info_file, data_info )


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

					time_adj = time_vol + conf[ "ana" ][ "hrf_corr_vol" ]

					run_info[ i_run, i_block, 0 ] = time_adj

	np.save( paths[ "mvpa" ][ "blk_data_info" ], run_info )


def blk_data_xtr( paths, conf ):
	"""Extracts blocked and z-scored data"""

	n_vol_per_block = conf[ "exp" ][ "block_len_s" ] / conf[ "acq" ][ "tr_s" ]

	blk_data_info = np.load( paths[ "mvpa" ][ "blk_data_info" ] )

	for hemi in [ "lh", "rh" ]:

		# data is runs x nodes x time
		data = np.load( "%s-%s.npy" % ( paths[ "mvpa" ][ "data" ], hemi ) )
		data_info = np.load( "%s-%s.npy" % ( paths[ "mvpa" ][ "data_info" ], hemi ) )

		( n_runs, n_nodes, n_vols ) = data.shape
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

				# this extraction flips the dimensions around; `blk_tc` is vol x node,
				# surprisingly
				blk_tc = data[ i_run, :, i_vol_range ]

				# just verify that is indeed the case
				assert( blk_tc.shape == ( len( i_vol_range ), n_nodes ) )

				# ... because we want to average over the vol dimension
				blk_mean = np.mean( blk_tc, axis = 0 )

				# average over block timepoints
				blk_data[ i_run, i_block, : ] = blk_mean

			# z-score all the blocks in this run, for each node
			blk_data[ i_run, :, : ] = scipy.stats.zscore( blk_data[ i_run, :, : ],
			                                              axis = 0
			                                            )

		blk_file = "%s-%s.npy" % ( paths[ "mvpa" ][ "blk_data" ], hemi )

		np.save( blk_file, blk_data )


def searchlight( paths, conf ):
	"""Runs a searchlight analysis"""

	# load the info for each block, which doesn't depend on hemisphere
	blk_data_info = np.load( paths[ "mvpa" ][ "blk_data_info" ] )

	svm_temp_dir = os.path.join( paths[ "mvpa" ][ "base_dir" ], "svm_temp" )

	train_info = [ [ 2, 3, 4, 5, 6, 7, 8, 9 ],
	               [ 4, 5, 6, 7, 8, 9, 0, 1 ],
	               [ 6, 7, 8, 9, 0, 1, 2, 3 ],
	               [ 8, 9, 0, 1, 2, 3, 4, 5 ],
	               [ 0, 1, 2, 3, 4, 5, 6, 7 ]
	             ]

	test_info = [ [ 0, 1 ],
	              [ 2, 3 ],
	              [ 4, 5 ],
	              [ 6, 7 ],
	              [ 8, 9 ]
	            ]

	# analysis proceeds separately for the two hemispheres
	for hemi in [ "lh", "rh" ]:

		# load the all-important data
		blk_data = np.load( "%s-%s.npy" % ( paths[ "mvpa" ][ "blk_data" ],
		                                    hemi
		                                  )
		                  )

		# ... and the data info; ( ( node index, roi ), node )
		data_info = np.load( "%s-%s.npy" % ( paths[ "mvpa" ][ "data_info" ],
		                                     hemi
		                                   )
		                   )

		( n_runs, n_blocks, n_nodes ) = blk_data.shape

		# first step is to transform the ROI dataset into a more readable form

		# we don't need (or want) the full version here, since we write out the
		# node indices
		roi_file = "%s_%s.niml.dset" % ( paths[ "rois" ][ "dset" ], hemi )
		seed_node_file = "%s_%s" % ( paths[ "mvpa" ][ "nodes" ], hemi )

		conv_cmd = [ "ConvertDset",
		             "-o_1D",  # output format, 1D
		             "-input", roi_file,
		             "-prepend_node_index_1D",
		             "-overwrite",
		             "-prefix", seed_node_file
		           ]

		fmri_tools.utils.run_cmd( conv_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# now load the resulting nodes file
		seed_nodes = np.loadtxt( "%s.1D.dset" % seed_node_file ).astype( "int" )

		# path to a text file containing the current node number
		curr_seed_file = os.path.join( paths[ "mvpa" ][ "base_dir" ],
		                               "curr_seed.1D"
		                             )

		# base path for the current disk ROI info
		curr_disk_file = os.path.join( paths[ "mvpa" ][ "base_dir" ],
		                               "curr_disk"
		                             )

		pbar = progressbar.ProgressBar( widgets = [ progressbar.Percentage(),
		                                            progressbar.Bar()
		                                          ],
		                                maxval = seed_nodes.shape[ 0 ]
		                              ).start()

		acc = np.empty( ( n_nodes ) )
		acc.fill( np.NAN )

		# iterate through each seed node
		for ( i_seed_node, seed_node ) in enumerate( seed_nodes[ :, 0 ] ):

			pbar.update( i_seed_node )

			# write the current node to a file
			np.savetxt( curr_seed_file, [ seed_node ] )

			grow_cmd = [ "ROIgrow",
			             "-i", "%s-%s.asc" % ( paths[ "reg" ][ "flat" ], hemi ),
			             "-roi_nodes", curr_seed_file,
			             "-lim", "%d" % conf[ "ana" ][ "slight_r" ],
			             "-prefix", curr_disk_file,
			             "-overwrite"
			           ]

			fmri_tools.utils.run_cmd( grow_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

			# load the nodes corresponding to this seed node
			disk_nodes = np.loadtxt( "%s.1D" % curr_disk_file )

			if disk_nodes.size == 1:
				disk_nodes = np.array( [ disk_nodes ] )

			disk_data = np.empty( ( n_runs, n_blocks, len( disk_nodes ) ) )
			disk_data.fill( np.NAN )

			i_disk_nodes = []

			for disk_node in disk_nodes:

				i = np.where( data_info[ 0, : ] == disk_node )[ 0 ]

				if i.size == 1:

					i_disk_node = i[ 0 ]

					if i_disk_node < n_nodes:
						i_disk_nodes.append( i_disk_node )

			i_disk_nodes = np.array( i_disk_nodes )

			disk_data = blk_data[ :, :, i_disk_nodes ]

			assert( disk_data.shape == ( n_runs, n_blocks, len( i_disk_nodes ) ) )

			# now we have our data set, can run the classification
			acc[ i_seed_node ] = _svm_classify( disk_data,
			                                    blk_data_info,
			                                    train_info,
			                                    test_info,
			                                    svm_temp_dir,
			                                    paths[ "summ" ][ "log_file" ]
			                                  )

		pbar.finish()

		np.savetxt( "%s-%s.txt" % ( paths[ "mvpa" ][ "acc" ], hemi ), acc )


def write_searchlight( paths, conf ):
	"""Write the searchlight analysis to surface datasets"""

	start_dir = os.getcwd()

	os.chdir( paths[ "mvpa" ][ "base_dir" ] )

	for hemi in [ "lh", "rh" ]:

		acc_file = "%s-%s.txt" % ( paths[ "mvpa" ][ "acc" ], hemi )

		seed_node_file = "%s_%s.1D.dset" % ( paths[ "mvpa" ][ "nodes" ], hemi )

		out_file = "%s-%s.niml.dset" % ( paths[ "mvpa" ][ "acc" ], hemi )

		conv_cmd = [ "ConvertDset",
		             "-o_niml",
		             "-input", acc_file,
		             "-node_index_1D", "%s[0]" % seed_node_file,
		             "-i_1D",
		             "-prefix", out_file
		           ]

		fmri_tools.utils.run_cmd( conv_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = log_path
		                        )

	os.chdir( start_dir )


def _svm_classify( blk_data,
                   blk_data_info,
                   train_info,
                   test_info,
                   temp_dir,
                   log_path
                 ):
	"""Helper function to run SVM classification"""

	( n_runs, n_blocks, n_nodes ) = blk_data.shape

	svm_runs = zip( train_info, test_info )

	n_folds = len( svm_runs )

	acc = np.empty( ( n_folds ) )
	acc.fill( np.NAN )

	for ( i_fold, ( train_runs, test_runs ) ) in enumerate( svm_runs ):

		train_data = blk_data[ train_runs, :, : ]

		( n_train_runs, n_train_blocks, n_train_nodes ) = train_data.shape

		test_data = blk_data[ test_runs, :, : ]

		( n_test_runs, n_test_blocks, n_test_nodes ) = test_data.shape

		train_path = os.path.join( temp_dir, "train.txt" )

		train_file = open( train_path, "w" )

		for i_run in xrange( n_train_runs ):
			for i_block in xrange( n_train_blocks ):

				ex_cond = blk_data_info[ train_runs[ i_run ], i_block, 1 ]

				ex_str = "%+d" % ex_cond

				for i_node in xrange( n_train_nodes ):

					ex_str += " %d:%.12f" % ( i_node + 1,
					                          train_data[ i_run, i_block, i_node ]
					                        )

				train_file.write( "%s\n" % ex_str )

		train_file.close()

		model_path = os.path.join( temp_dir, "model.txt" )

		train_cmd = [ "svm_learn",
		              train_path,
		              model_path
		            ]

		fmri_tools.utils.run_cmd( train_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = log_path
		                        )

		test_path = os.path.join( temp_dir, "test.txt" )

		test_file = open( test_path, "w" )

		true_conds = []

		for i_run in xrange( n_test_runs ):
			for i_block in xrange( n_test_blocks ):

				ex_cond = blk_data_info[ test_runs[ i_run ], i_block, 1 ]

				true_conds.append( ex_cond )

				ex_str = "%+d" % ex_cond

				for i_node in xrange( n_test_nodes ):

					ex_str += " %d:%.12f" % ( i_node + 1,
					                          test_data[ i_run, i_block, i_node ]
					                        )

				test_file.write( "%s\n" % ex_str )

		true_conds = np.array( true_conds )

		test_file.close()

		pred_path = os.path.join( temp_dir, "pred.txt" )

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

		acc[ i_fold ] = np.sum( true_conds == np.sign( pred ) ) / len( pred ) * 100

	return np.mean( acc )


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
