"""
Set of routines to analyse single-subject fMRI data for the natural scenes
apertures fMRI experiment.
"""

from __future__ import division

import os, os.path
import logging

import numpy as np
import scipy.stats
import svmlight
import progressbar

import fmri_tools.utils

import ns_aperture.fmri.exp, ns_aperture.fmri.loc


def design_prep( conf, paths ):
	"""Prepares the designs for GLM analysis"""

	seq_info = ns_aperture.fmri.exp.get_seq_ind()

	# open the condition file to write
	cond_file = open( paths.ana.stim_times.full( ".txt" ), "w" )

	for run_num in conf.subj.exp_runs:

		run_seq = np.load( paths.exp_log.seq.full( "_{n:d}.npy".format( n = run_num ) ) )

		( n_evt, _ ) = run_seq.shape

		for i_evt in xrange( n_evt ):

			curr_block_num = run_seq[ i_evt, seq_info[ "block_num" ] ]
			prev_block_num = run_seq[ i_evt - 1, seq_info[ "block_num" ] ]

			is_transition = ( curr_block_num != prev_block_num )

			if is_transition:

				start_time_s = run_seq[ i_evt, seq_info[ "time_s" ] ]

				cond = int( run_seq[ i_evt, seq_info[ "block_type" ] ] )

				# don't want to model non-coherent blocks
				if ( cond == 0 ):
					cond_file.write( "{x:.03f} ".format( x = start_time_s ) )

					assert ( run_seq[ i_evt, seq_info[ "img_id_L" ] ] ==
					         run_seq[ i_evt, seq_info[ "img_id_R" ] ]
					       )

		cond_file.write( "\n" )

	cond_file.close()


def glm( conf, paths, std_surf = True ):
	"""Experiment GLM"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running GLM..." )

	n_vols_per_block = int( conf.exp.block_len_s / conf.acq.tr_s )
	n_vols = int( conf.exp.run_len_s / conf.acq.tr_s )
	cens_str = "*:0-{n:d}".format( n = n_vols_per_block - 1 )
	cens_str += " *:{p:d}-{n:d}".format( p = n_vols - n_vols_per_block,
	                                     n = n_vols - 1
	                                   )

	start_dir = os.getcwd()

	os.chdir( paths.ana.base.full() )

	exp_surfs = [ paths.func.surfs[ i - 1 ] for i in conf.subj.exp_runs ]

	for hemi in [ "lh", "rh" ]:

		if std_surf:
			hemi_ext = "-std_{h:s}".format( h = hemi )
		else:
			hemi_ext = "_{h:s}".format( h = hemi )

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		if std_surf:
			surf_paths = [ surf_path.full( "-smooth{h:s}.niml.dset".format( h = hemi_ext ) )
			               for surf_path in exp_surfs
			             ]
		else:
			surf_paths = [ surf_path.full( "{h:s}.niml.dset".format( h = hemi_ext ) )
			               for surf_path in exp_surfs
			             ]

		glm_cmd.extend( surf_paths )

		glm_cmd.extend( [ "-force_TR", "{tr:.3f}".format( tr = conf.acq.tr_s ),
		                  "-polort", conf.ana.poly_ord,
		                  "-local_times",
		                  "-CENSORTR", cens_str,
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-overwrite",
		                  "-x1D_stop",  # want to use REML, so don't bother running
		                  "-num_stimts", "1",
		                  "-stim_label", "1", "coh",
		                  "-stim_times", "1",
		                                 paths.ana.stim_times.full( ".txt" ),
		                                 conf.ana.hrf_model,
		                  "-gltsym", "'SYM: +coh'",
		                  "-glt_label", "1", "coh_gt"
		                ]
		              )

		# run this first GLM
		fmri_tools.utils.run_cmd( " ".join( glm_cmd ) )

		# delete the annoying command file that 3dDeconvolve writes
		os.remove( "Decon.REML_cmd" )

		beta_file = paths.ana.beta.file( hemi_ext + ".niml.dset" )
		buck_file = paths.ana.glm.file( hemi_ext + ".niml.dset" )

		reml_cmd = [ "3dREMLfit",
		             "-matrix", "exp_design.xmat.1D",
		             "-Rbeta", beta_file,
		             "-tout",
		             "-Rbuck", buck_file,
		             "-overwrite",
		             "-input"
		           ]

		reml_cmd.append( "'" + " ".join( surf_paths ) + "'" )

		# run the proper GLM
		fmri_tools.utils.run_cmd( " ".join( reml_cmd ) )

		if std_surf:

			in_dsets = [ beta_file, buck_file ]
			out_dsets = [ paths.ana.beta.file( hemi_ext + "-full.niml.dset" ),
			              paths.ana.glm.file( hemi_ext + "-full.niml.dset" )
			            ]

			for ( in_dset, out_dset ) in zip( in_dsets, out_dsets ):
				# convert the beta and glm files to full
				fmri_tools.utils.sparse_to_full( in_dset = in_dset,
				                                 out_dset = out_dset,
				                                 pad_node = "ld141"
				                               )


	os.chdir( start_dir )


def coh_clust_summ( conf, paths ):
	"Summarise the significant cluster activation"

	group_conf = ns_aperture.config.get_conf()
	group_paths = ns_aperture.paths.get_group_paths( group_conf )

	os.chdir( paths.ana.base.full() )

	for hemi in [ "lh", "rh" ]:

		hemi_ext = "-std_{h:s}".format( h = hemi )

		beta_path = paths.ana.glm.file( hemi_ext + "-full.niml.dset'[3]'" )

		out_path = paths.ana.clust.file( hemi_ext + ".txt" )

		clust_path = group_paths.coh_clust.full( hemi_ext + "-full.niml.dset" )

		cmd = [ "3dmaskdump",
		        "-mask", clust_path,
		        "-noijk",
		        "-o", out_path,
		        beta_path
		      ]

		fmri_tools.utils.run_cmd( " ".join( cmd ) )



def loc_design_prep( conf, paths ):
	"""Prepares the designs for GLM analysis"""

	seq_info = ns_aperture.fmri.loc.get_seq_ind()

	# L /R
	n_cond = 2

	cond_files = [ open( "%s_%d.txt" % ( paths.loc.time_files,
	                                    cond_num
	                                  ),
	                     "w"
	                   )
	                   for cond_num in np.arange( 1, n_cond + 1 )
	             ]

	loc_ord = ( "AB", "BA" )

	for i_run in xrange( conf.subj.n_loc_runs ):

		run_times = []
		run_conds = []

		run_seq = ns_aperture.fmri.loc.get_seq( conf, loc_ord[ i_run ] )

		( n_evt, _ ) = run_seq.shape

		for i_evt in xrange( n_evt ):

			curr_block_num = run_seq[ i_evt, seq_info[ "block_num" ] ]
			prev_block_num = run_seq[ i_evt - 1, seq_info[ "block_num" ] ]

			is_transition = ( curr_block_num != prev_block_num )

			if is_transition:

				start_time_s = run_seq[ i_evt, seq_info[ "time_s" ] ]

				cond = int( run_seq[ i_evt, seq_info[ "block_type" ] ] )

				# 0 = blank
				if cond > 0:
					run_times.append( start_time_s )
					run_conds.append( cond - 1 )

		run_times = np.array( run_times )
		run_conds = np.array( run_conds )

		for i_cond in xrange( n_cond ):

			i_evt_cond = np.where( run_conds == i_cond )[ 0 ]

			if i_evt_cond.size == 0:
				cond_files[ i_cond ].write( "*" )
			else:
				_ = [ cond_files[ i_cond ].write( "%.5f\t" % evt_time )
				      for evt_time in run_times[ i_evt_cond ]
				    ]

			cond_files[ i_cond ].write( "\n" )

	_ = [ cond_file.close() for cond_file in cond_files ]


def loc_glm( conf, paths, std_surf = True ):
	"""Loclaiser GLM"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running localiser GLM..." )

	block_len_s = int( conf.exp.loc_run_len_s / conf.exp.loc_n_blocks )
	n_vols_per_block = int( block_len_s / conf.acq.tr_s )
	n_vols = int( conf.exp.loc_run_full_len_s / conf.acq.tr_s )
	n_pre_vols = int( conf.exp.loc_pre_len_s / conf.acq.tr_s )
	cens_str = "*:0-{n:d}".format( n = n_pre_vols - 1 )

	n_cond = 2

	hrf_model = "SPMG1({n:d})".format( n = block_len_s )

	stim_files = [ "%s_%d.txt" % ( paths.loc.time_files,
	                              cond_num
	                            )
	               for cond_num in np.arange( 1, n_cond + 1 )
	             ]

	stim_labels = ( "L", "R" )


	start_dir = os.getcwd()

	os.chdir( paths.loc.base.full() )

	exp_surfs = [ paths.func.surfs[ i - 1 ] for i in conf.subj.loc_runs ]

	for hemi in [ "lh", "rh" ]:

		if std_surf:
			hemi_ext = "-std_{h:s}".format( h = hemi )
		else:
			hemi_ext = "_{h:s}".format( h = hemi )

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		if std_surf:
			surf_paths = [ surf_path.full( "-smooth{h:s}.niml.dset".format( h = hemi_ext ) )
			               for surf_path in exp_surfs
			             ]
		else:
			surf_paths = [ surf_path.full( "{h:s}.niml.dset".format( h = hemi_ext ) )
			               for surf_path in exp_surfs
			             ]

		glm_cmd.extend( surf_paths )

		glm_cmd.extend( [ "-force_TR", "{tr:.3f}".format( tr = conf.acq.tr_s ),
		                  "-polort", conf.ana.poly_ord,
		                  "-local_times",
		                  "-CENSORTR", cens_str,
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-overwrite",
		                  "-x1D_stop",  # want to use REML, so don't bother running
		                  "-num_stimts", "2"
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

		glm_cmd.extend( [ "-gltsym", "'SYM: +L +R'", "-glt_label", "1", "both",
		                  "-gltsym", "'SYM: +L -R'", "-glt_label", "2", "lGTr"
		                ]
		              )

		# run this first GLM
		fmri_tools.utils.run_cmd( " ".join( glm_cmd ) )

		# delete the annoying command file that 3dDeconvolve writes
		os.remove( "Decon.REML_cmd" )

		beta_file = paths.loc.beta.file( hemi_ext + ".niml.dset" )
		buck_file = paths.loc.glm.file( hemi_ext + ".niml.dset" )

		reml_cmd = [ "3dREMLfit",
		             "-matrix", "exp_design.xmat.1D",
		             "-Rbeta", beta_file,
		             "-tout",
		             "-Rbuck", buck_file,
		             "-overwrite",
		             "-input"
		           ]

		reml_cmd.append( "'" + " ".join( surf_paths ) + "'" )

		# run the proper GLM
		fmri_tools.utils.run_cmd( " ".join( reml_cmd ) )

		if std_surf:

			in_dsets = [ beta_file, buck_file ]
			out_dsets = [ paths.loc.beta.file( hemi_ext + "-full.niml.dset" ),
			              paths.loc.glm.file( hemi_ext + "-full.niml.dset" )
			            ]

			for ( in_dset, out_dset ) in zip( in_dsets, out_dsets ):
				# convert the beta and glm files to full
				fmri_tools.utils.sparse_to_full( in_dset = in_dset,
				                                 out_dset = out_dset,
				                                 pad_node = "ld141"
				                               )


	os.chdir( start_dir )


def mvpa_prep( conf, paths ):
	"Prepare the structures for MVPA analysis"

	# Needs:
	#  1. List of node indices for each hemisphere
	#  2. Block summary data for each node, run, and block for each hemisphere
	#  3. Condition info for each run and block

	grp_paths = ns_aperture.paths.get_group_paths( conf )

	# starts with a coherent block, but the first block is discarded so actually
	# begins with a non-coherent block
	cond_info = np.repeat( np.tile( [ -1, +1 ],
	                                int( conf.ana.mvpa_blocks / 2 )
	                              )[ np.newaxis, : ],
	                       conf.subj.n_exp_runs,
	                       axis = 0
	                     )

	# flip alternate runs
	cond_info[ 1::2, : ] = cond_info[ 1::2, ::-1 ]

	cond_path = paths.mvpa.cond_info.full( ".txt" )

	np.savetxt( cond_path, cond_info, "%d" )

	for hemi in [ "lh", "rh" ]:

		data = np.empty( ( conf.ana.n_common_nodes[ hemi ],
		                   len( conf.subj.exp_runs ),
		                   conf.ana.mvpa_blocks
		                 )
		               )
		data.fill( np.NAN )

		mask_path = grp_paths.surf_mask.full( "-std_" + hemi + "-full.niml.dset" )

		for ( i_run, run_num ) in enumerate( conf.subj.exp_runs ):

			# extract the data
			surf = paths.func.filts[ run_num - 1 ]

			os.chdir( surf.dir() )

			surf_file = surf.file( "-std_{h:s}-full.niml.dset".format( h = hemi ) )

			out_file = surf.file( "-std_{h:s}-full.txt".format( h = hemi ) )

			if os.path.exists( out_file ):
				os.remove( out_file )

			cmd = [ "3dmaskdump",
			        "-mask", mask_path,
			        "-index",
			        "-noijk",
			        "-nozero",
			        "-o", out_file,
			        surf_file
			      ]

			fmri_tools.utils.run_cmd( " ".join( cmd ) )

			run_data = np.loadtxt( out_file )

			# only need to save the nodes once
			if i_run == 0:

				nodes = run_data[ :, 0 ]
				node_path = paths.mvpa.nodes.full( "_" + hemi + ".txt" )
				np.savetxt( node_path, nodes, "%d" )

			# just check that the nodes are all the same across runs
			assert np.all( nodes == run_data[ :, 0 ] )

			# chop out the nodes
			run_data = run_data[ :, 1: ]

			for i_block in xrange( conf.ana.mvpa_blocks ):

				i_start = conf.ana.block_len_vol * i_block + conf.ana.i_start_hrf_corr

				i_block_vols = np.arange( i_start, i_start + conf.ana.block_len_vol )

				data[ :, i_run, i_block ] = np.mean( run_data[ :, i_block_vols ],
				                                     axis = 1
				                                   )


		# z-score
		data = scipy.stats.zscore( data, axis = 2 )

		data_path = paths.mvpa.data.full( "_" + hemi + ".npy" )

		np.save( data_path, data )


def mvpa( conf, paths ):
	"Runs the MVPA analysis"

	group_conf = ns_aperture.config.get_conf()
	group_paths = ns_aperture.paths.get_group_paths( group_conf )

	cond_info = np.loadtxt( paths.mvpa.cond_info.full( ".txt" ), np.int )

	for hemi in [ "lh", "rh" ]:

		# nodes x runs x blocks
		data = np.load( paths.mvpa.data.full( "_" + hemi + ".npy" ) )

		seed_nodes = np.loadtxt( paths.mvpa.nodes.full( "_" + hemi + ".txt" ),
		                         np.int
		                       )

		acc = np.empty( ( data.shape[ 0 ] ) )
		acc.fill( np.NAN )

		pbar = progressbar.ProgressBar( widgets = [ progressbar.Percentage(),
		                                            progressbar.Bar()
		                                          ],
		                                maxval = seed_nodes.shape[ 0 ]
		                              ).start()

		with open( group_paths.sl_info.full( "_" + hemi + ".txt" ), "r" ) as sl_info:

			for ( i_seed, ( seed_node, node_line ) ) in enumerate( zip( seed_nodes, sl_info.readlines() ) ):

				pbar.update( i_seed )

				nodes = [ int( x ) for x in node_line.splitlines()[ 0 ].split( "\t" ) ]

				assert seed_node in nodes

				i_nodes = []

				for sl_node in nodes:
					i = np.where( seed_nodes == sl_node )[ 0 ]
					if len( i ) > 0:
						assert len( i ) == 1
						i_nodes.append( i[ 0 ] )

				if len( i_nodes ) > 0:
					acc[ i_seed ] = _classify( data[ i_nodes, ... ], cond_info )

		# save acc
		acc_path = paths.mvpa.acc.full( "_" + hemi + ".txt" )
		np.savetxt( acc_path, acc )

		pbar.finish()

		os.chdir( paths.mvpa.base.full() )

		# convert to full niml
		cmd = [ "ConvertDset",
		        "-i_1D", "-input", acc_path,
		        "-node_index_1D", paths.mvpa.nodes.full( "_" + hemi + ".txt" ),
		        "-o_niml", "-prefix", paths.mvpa.acc.full( "_" + hemi + ".niml.dset" ),
		        "-pad_to_node", "ld141",
		        "-overwrite"
		      ]

		fmri_tools.utils.run_cmd( " ".join( cmd ) )


def _classify( data, cond_info ):
	"Runs a single classification"

	( n_runs, n_blocks ) = cond_info.shape

	acc = np.empty( ( n_runs ) )
	acc.fill( np.NAN )

	for i_test_run in xrange( n_runs ):

		# exclude the test run from the training set
		i_train = np.setdiff1d( np.arange( n_runs ), [ i_test_run ] )

		train_data = _format_data( data[ :, i_train, : ],
		                           cond_info[ i_train, : ]
		                         )

		model = svmlight.learn( train_data,
		                        type = "classification",
		                        kernel = "linear"
		                      )

		test_data = _format_data( data[ :, [ i_test_run ], : ],
		                          cond_info[ [ i_test_run ], : ]
		                        )

		pred = svmlight.classify( model, test_data )

		f_acc = ( float( ( np.sign( pred ) == cond_info[ i_test_run, : ] ).sum() ) /
		          len( pred ) * 100.0
		        )

		acc[ i_test_run ] = f_acc

	assert np.sum( np.isnan( acc ) ) == 0

	return np.mean( acc )


def _format_data( data, cond_info ):

	f_data = []

	for i_run in xrange( data.shape[ 1 ] ):
		for i_block in xrange( data.shape[ 2 ] ):

			eg_data = []

			for i_node in xrange( data.shape[ 0 ] ):

				eg_data.extend( [ ( i_node + 1, data[ i_node, i_run, i_block ] ) ] )

			f_data.append( ( cond_info[ i_run, i_block ], eg_data ) )

	return f_data


def write_searchlight( conf, paths, node_num, hemi, base_path ):
	"Debugging script to make sure searchlights look as expected"

	group_conf = ns_aperture.config.get_conf()
	group_paths = ns_aperture.paths.get_group_paths( group_conf )

	( base_dir, base_file ) = os.path.split( base_path )

	os.chdir( base_dir )

	node_found = False

	node_list = np.loadtxt( paths.mvpa.nodes.full( "_" + hemi + ".txt" ) )

	with open( group_paths.sl_info.full( "_" + hemi + ".txt" ), "r" ) as sl_info:

		for ( seed_node, node_line ) in zip( node_list, sl_info.readlines() ):

			if seed_node != node_num:
				continue

			nodes = [ int( x ) for x in node_line.splitlines()[ 0 ].split( "\t" ) ]

			node_found = True

			nodes_path = base_file + "-nodes.txt"

			np.savetxt( nodes_path, nodes, fmt = "%d" )

			data_path = base_file + "-data.txt"

			np.savetxt( data_path, np.ones( len( nodes ) ), fmt = "%d" )

			cmd = [ "ConvertDset",
			        "-i_1D", "-input", data_path,
			        "-node_index_1D", nodes_path,
			        "-o_niml", "-prefix", base_file + "-searchlight.niml.dset",
			        "-overwrite"
			      ]

			fmri_tools.utils.run_cmd( " ".join( cmd ) )


	if not node_found:
		print "Error: node not found"


def task_xtr( conf, paths ):
	"""Extracts the details on task performance"""

	task = []

	# different items are different 'hand's
	keys = [ { "r" : 0,
	           "g" : 1,
	           "y" : 2,
	           "b" : 3
	         },
	         { "r" : 3,
	           "g" : 2,
	           "y" : 1,
	           "b" : 0
	         }
	       ]

	# order was ordered for subject s1000, and the order was not stored in the
	# task file
	s1000_hand_order = [ 1, 0, 0, 1, 1, 0, 0, 1, 1, 0 ]

	for i_run in xrange( conf.subj.n_exp_runs ):

		task_file = paths.exp_log.task.full( "_" + str( i_run + 1 ) + ".npy" )
		run_tasks = np.load( task_file )

		seq_file = paths.exp_log.seq.full( "_" + str( i_run + 1 ) + ".npy" )
		run_seq = np.load( seq_file )

		# potential times that a task could have been cued
		evt_times = run_seq[ :, 0 ] + conf.task.evt_on_s

		for curr_task in run_tasks:

			if conf.subj.subj_id == "s1000":
				( resp, time_s ) = curr_task
				hand_flag = s1000_hand_order[ i_run ]
			else:
				( resp, time_s, hand_flag ) = curr_task

			# s1032 had it backwards
			if conf.subj.subj_id == "s1032":
				hand_flag = np.logical_not( hand_flag )

			# find the last event where the response time is greater than the
			# presentation onset time
			i_evt = np.where( time_s > evt_times )[ 0 ][ -1 ]

			( cond, img_id_L, img_id_R ) = run_seq[ i_evt, [ 2, 5, 6 ] ]

			resp_time = time_s - run_seq[ i_evt, 0 ]

			task.append( [ i_run,
			               cond,
			               keys[ hand_flag ][ resp ],
			               hand_flag,
			               img_id_L,
			               img_id_R,
			               resp_time
			             ]
			           )

	np.save( paths.ana.task.full( ".npy" ), task )

