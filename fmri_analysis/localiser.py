"""
Set of routines to analyse single-subject localiser fMRI data for the natural
scenes apertures fMRI experiment.
"""

from __future__ import division

import os, os.path

import numpy as np
import scipy.stats

import fmri_tools.utils

import ns_aperture.fmri.loc


def loc_design( paths, conf ):
	"""Prepares the designs for GLM analysis"""

	seq_info = ns_aperture.fmri.loc.get_seq_ind()

	# L /R
	n_cond = 2

	cond_files = [ open( "%s%d.txt" % ( paths[ "loc" ][ "time_files" ],
	                                    cond_num
	                                  ),
	                     "w"
	                   )
	               for cond_num in np.arange( 1, n_cond + 1 )
	             ]

	loc_ord = ( "AB", "BA" )

	for i_run in xrange( conf[ "subj" ][ "n_loc_runs" ] ):

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


def loc_glm( paths, conf ):
	"""Localiser GLM"""

	start_dir = os.getcwd()

	os.chdir( paths[ "loc" ][ "base_dir" ] )

	n_cond = 2

	hrf_model = conf[ "ana" ][ "hrf_model" ]

	stim_files = [ "%s%d.txt" % ( paths[ "loc" ][ "time_files" ],
	                              cond_num
	                            )
	               for cond_num in np.arange( 1, n_cond + 1 )
	             ]

	stim_labels = ( "L", "R" )

	for hemi in [ "lh", "rh" ]:

		glm_cmd = [ "3dDeconvolve",
		            "-input"
		          ]

		i_surf_files = np.array( conf[ "subj" ][ "loc_runs" ] ).astype( "int" ) - 1

		surf_files = [ paths[ "func" ][ "surf_files" ][ i_surf ]
		               for i_surf in i_surf_files
		             ]

		glm_cmd.extend( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                  for surf_file in surf_files
		                ]
		              )

		glm_cmd.extend( [ "-force_TR", "%.3f" % conf[ "acq" ][ "tr_s" ],
		                  "-polort", "%d" % ( conf[ "ana" ][ "poly_ord" ] - 1 ),
		                  "-local_times",
		                  "-xjpeg", "exp_design.png",
		                  "-x1D", "exp_design",
		                  "-overwrite",
		                  "-x1D_stop",  # want to use REML, so don't bother running
		                  "-num_stimts", "%d" % n_cond
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

		# run this first GLM
		fmri_tools.utils.run_cmd( glm_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# delete the annoying command file that 3dDeconvolve writes
		os.remove( "Decon.REML_cmd" )

		glm_file = "%s_%s" % ( paths[ "loc" ][ "glm" ], hemi )

		reml_cmd = [ "3dREMLfit",
		             "-matrix", "exp_design.xmat.1D",
		             "-tout",
		             "-Rbuck", "%s_reml.niml.dset" % glm_file,
		             "-overwrite",
		             "-input"
		           ]

		reml_cmd.append( " ".join( [ "%s_%s.niml.dset" % ( surf_file, hemi )
		                             for surf_file in surf_files
		                           ]
		                         )
		               )

		# run the proper GLM
		fmri_tools.utils.run_cmd( reml_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# convert to full
		full_glm_file = "%s_%s_reml-full.niml.dset" % ( paths[ "loc" ][ "glm" ], hemi )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		fmri_tools.utils.sparse_to_full( "%s_reml.niml.dset" % glm_file,
		                                 full_glm_file,
		                                 pad_node = pad_node,
		                                 log_path = paths[ "summ" ][ "log_file" ],
		                                 overwrite = True
		                               )

	os.chdir( start_dir )


def loc_parcellate( paths, conf ):
	"""Extract t statistics from rois"""

	# LVF, RVF
	loc_stat_bricks = "2,4"

	os.chdir( paths[ "loc" ][ "base_dir" ] )

	# check by running:
	#  3dAttribute BRICK_STATSYM file
	df = 284

	# get the critical t statistic value
	crit_t = scipy.stats.t.isf( conf[ "ana" ][ "loc_p" ], df )

	for hemi in [ "lh", "rh" ]:

		loc_stat_file = "%s_%s_reml-full.niml.dset[%s]" % ( paths[ "loc" ][ "glm" ],
		                                                    hemi,
		                                                    loc_stat_bricks
		                                                  )

		roi_file = "%s_%s-full.niml.dset" % ( paths[ "rois" ][ "dset" ], hemi )

		for ( roi_name, roi_val ) in conf[ "ana" ][ "rois" ]:

			roi_stat_file = "%s_%s_%s.txt" % ( paths[ "loc" ][ "roi_stat" ],
			                                   roi_name,
			                                   hemi
			                                 )

			xtr_cmd = [ "3dmaskdump",
			            "-mask", roi_file,
			            "-mrange", roi_val, roi_val,
			            "-o", roi_stat_file,
			            "-index",  # save the node index
			            "-noijk",  # j and k are just zeros
			            loc_stat_file
			          ]

			# 3dmaskdump won't overwrite, so need to manually remove any prior file
			if os.path.exists( roi_stat_file ):
				os.remove( roi_stat_file )

			fmri_tools.utils.run_cmd( xtr_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

			roi_stat = np.loadtxt( roi_stat_file )


			# add a column of zeros that will hold the parcel number
			roi_stat = np.hstack( ( roi_stat,
			                        np.zeros( ( roi_stat.shape[ 0 ], 1 ) )
			                      )
			                    )

			# index 0 is the node index
			lvf = [ roi_stat[ :, 1 ] > crit_t,
			        np.abs( roi_stat[ :, 1 ] ) <= crit_t,
			        roi_stat[ :, 1 ] < -crit_t
			      ]

			rvf = [ roi_stat[ :, 2 ] > crit_t,
			        np.abs( roi_stat[ :, 2 ] ) <= crit_t,
			        roi_stat[ :, 2 ] < -crit_t
			      ]

			i_parcel = 1

			for i_lvf in xrange( 3 ):
				for i_rvf in xrange( 3 ):

					in_parcel = np.logical_and( lvf[ i_lvf ], rvf[ i_rvf ] )

					assert( np.all( roi_stat[ in_parcel, 3 ] == 0 ) )

					roi_stat[ in_parcel, 3 ] = i_parcel

					i_parcel += 1

			assert( np.logical_not( np.any( roi_stat[ :, 3 ] == 0 ) ) )

			np.savetxt( roi_stat_file, roi_stat )

			roi_parc_file = "%s_%s_%s-full.niml.dset" % ( paths[ "loc" ][ "roi_parc" ],
			                                              roi_name,
			                                              hemi
			                                            )

			# now to convert it to a full niml dataset
			conv_cmd = [ "ConvertDset",
			             "-i_1D",
			             "-input", "%s[3]" % roi_stat_file,
			             "-o_niml",
			             "-prefix", roi_parc_file,
			             "-overwrite",
			             "-node_index_1D", "%s[0]" % roi_stat_file,
			             "-pad_to_node", "%d" % conf[ "subj" ][ "node_k" ][ hemi ]
			           ]

			fmri_tools.utils.run_cmd( conv_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )

	for ( roi_name, _ ) in conf[ "ana" ][ "rois" ]:

		parcel_file = "%s_%s.txt" % ( paths[ "loc" ][ "parc" ], roi_name )

		parcel_count = np.zeros( ( 3, 3 ) )

		for hemi in [ "lh", "rh" ]:

			roi_stat_file = "%s_%s_%s.txt" % ( paths[ "loc" ][ "roi_stat" ],
			                                   roi_name,
			                                   hemi
			                                 )

			roi_stat = np.loadtxt( roi_stat_file )

			i_parcel = 1

			for i_lvf in xrange( 3 ):
				for i_rvf in xrange( 3 ):

					parcel_count[ i_lvf, i_rvf ] += np.sum( roi_stat[ :, 3 ] == i_parcel )

					i_parcel += 1

		np.savetxt( parcel_file, parcel_count, "%d" )
