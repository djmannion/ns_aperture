"""
Set of routines to analyse single-subject fMRI data for the natural scenes
apertures fMRI experiment.
"""

from __future__ import division

import os, os.path
import logging

import numpy as np

import fmri_tools.utils

import ns_aperture.fmri.exp


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

	n_vols = int( conf.exp.run_len_s / conf.acq.tr_s )

	cens = np.ones( n_vols )

	n_vols_per_block = int( conf.exp.block_len_s / conf.acq.tr_s )

	cens[ :n_vols_per_block ] = 0
	cens[ -n_vols_per_block: ] = 0

	assert( np.sum( cens == 0 ) == ( n_vols_per_block * 2 ) )

	cens = np.tile( cens[ :, np.newaxis ], conf.subj.n_exp_runs ).T

	np.savetxt( paths.ana.cens.full( ".txt" ), cens, fmt = "%d" )


def glm( conf, paths, std_surf = False ):
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



def beta_to_psc( paths, conf ):
	"""Convert the GLM beta weights into units of percent signal change"""

	# this is the index in the beta file for the data we want to convert
	beta_brick = "0"

	start_dir = os.getcwd()

	os.chdir( paths[ "ana" ][ "base_dir" ] )

	for hemi in [ "lh", "rh" ]:

		# dataset holding the beta weights
		beta_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "beta" ], hemi )

		# design matrix file
		mat_file = os.path.join( paths[ "ana" ][ "base_dir" ],
		                         "exp_design.xmat.1D"
		                       )

		# baseline timecourse dataset, to write
		bltc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "bltc" ], hemi )

		# generate an average baseline timecourse
		bl_cmd = [ "3dSynthesize",
		           "-cbucket", beta_file,
		           "-matrix", mat_file,
		           "-cenfill", "none",  # this is important
		           "-select", "poly",  # only use the polynomials
		           "-prefix", bltc_file,
		           "-overwrite"
		         ]

		fmri_tools.utils.run_cmd( bl_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# baseline (point-estimate) dataset, to write
		bl_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "bl" ], hemi )

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
		psc_file = "%s_%s.niml.dset" % ( paths[ "ana" ][ "psc" ], hemi )

		# the input beta file, with sub-brick selector
		beta_sel = "%s[%s]" % ( beta_file, beta_brick )

		# check that the label is as expected
		beta_label = fmri_tools.utils.get_dset_label( beta_sel )
		assert( beta_label == [ "coh#0" ] )

		# compute psc
		# from http://afni.nimh.nih.gov/sscc/gangc/TempNorm.html
		psc_cmd = [ "3dcalc",
		            "-fscale",
		            "-a", bl_file,
		            "-b", beta_sel,
		            "-expr", "100 * b/a * step (1- abs(b/a))",
		            "-prefix", psc_file,
		            "-overwrite"
		          ]

		fmri_tools.utils.run_cmd( psc_cmd,
		                          env = fmri_tools.utils.get_env(),
		                          log_path = paths[ "summ" ][ "log_file" ]
		                        )

		# convert to full
		full_psc_file = "%s_%s-full" % ( paths[ "ana" ][ "psc" ], hemi )

		pad_node = "%d" % conf[ "subj" ][ "node_k" ][ hemi ]

		fmri_tools.utils.sparse_to_full( psc_file,
		                                 full_psc_file,
		                                 pad_node = pad_node,
		                                 log_path = paths[ "summ" ][ "log_file" ],
		                                 overwrite = True
		                               )

	os.chdir( start_dir )


def roi_xtr( paths, conf ):
	"""Extract PSC and statistics data from ROIs"""

	os.chdir( paths[ "rois" ][ "base_dir" ] )

	for ( roi_name, _ ) in conf[ "ana" ][ "rois" ]:

		for hemi in [ "lh", "rh" ]:

			# the *full* localiser mask file
			loc_mask_file = "%s_%s_%s-full.niml.dset" % ( paths[ "loc" ][ "roi_parc" ],
			                                              roi_name,
			                                              hemi
			                                            )

			roi_psc_file = "%s_%s_%s.txt" % ( paths[ "rois" ][ "psc" ],
			                                  roi_name,
			                                  hemi
			                                )

			# 3dmaskdump won't overwrite, so need to manually remove any prior file
			if os.path.exists( roi_psc_file ):
				os.remove( roi_psc_file )

			# our input dataset
			data_file = "%s_%s-full.niml.dset" % ( paths[ "ana" ][ "psc" ],
			                                       hemi
			                                     )

			# use the ROI file to mask the input dataset
			xtr_cmd = [ "3dmaskdump",
			            "-mask", loc_mask_file,
			            "-noijk",
			            "-o", roi_psc_file,
			            data_file,
			            loc_mask_file
			          ]

			fmri_tools.utils.run_cmd( xtr_cmd,
			                          env = fmri_tools.utils.get_env(),
			                          log_path = paths[ "summ" ][ "log_file" ]
			                        )


def roi_mean( paths, conf ):
	"""Average nodes within each parcel of each ROI"""

	os.chdir( paths[ "rois" ][ "base_dir" ] )

	for ( roi_name, _ ) in conf[ "ana" ][ "rois" ]:

		psc = np.vstack( [ np.loadtxt( "%s_%s_%s.txt" % ( paths[ "rois" ][ "psc" ],
		                                                  roi_name,
		                                                  hemi
		                                                )
		                             )
		                   for hemi in [ "lh", "rh" ]
		                 ]
		               )

		parc_mean = np.empty( len( conf[ "ana" ][ "i_parc" ] ) )
		parc_mean.fill( np.NAN )

		for ( i_parc, parc_ind ) in enumerate( conf[ "ana" ][ "i_parc" ] ):

			# psc array indices as those with a parcel number equal to any in the
			# list
			i_parc_psc = np.any( [ psc[ :, 1 ] == k
			                       for k in parc_ind
			                     ],
			                     axis = 0
			                   )

			parc_mean[ i_parc ] = np.mean( psc[ i_parc_psc, 0 ] )

		np.savetxt( "%s_%s.txt" % ( paths[ "rois" ][ "parc_psc" ],
		                            roi_name
		                          ),
		            parc_mean
		          )


def group_agg( grp_paths, conf ):
	"""Aggregate the univariate analyses across subjects"""

	n_subj = len( conf[ "all_subj" ] )

	n_rois = len( conf[ "ana" ][ "rois" ] ) * len( conf[ "ana" ][ "i_parc" ] )

	uni_data = np.empty( n_subj, n_rois )
	uni_data.fill( np.NAN )

	for ( i_subj, subj_id ) in conf[ "all_subj" ]:

		subj_conf = ns_aperture.config.get_conf( subj_id )
		subj_paths = ns_aperture.fmri_analysis.paths.get_subj_paths( subj_conf )

		i_roi = 0

		for ( roi_name, _ ) in conf[ "ana" ][ "rois" ]:

			roi_psc = np.loadtxt( paths[ "rois" ][ "parc_psc" ], roi_name )

			uni_data[ i_subj, i_roi:( i_roi + 2 ) ] = roi_psc

			i_roi += 2

	np.savetxt( grp_paths[ "uni" ], uni_data )
