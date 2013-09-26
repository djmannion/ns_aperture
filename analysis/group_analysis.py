
import os.path
import string

import numpy as np
import progressbar

import fmri_tools.utils, fmri_tools.preproc

import ns_aperture.config, ns_aperture.paths


def group_mask( conf, paths ):

	subj_ids = conf.all_subj.subj.keys()

	for hemi in [ "lh", "rh" ]:

		hemi_ext = "-std_" + hemi + "-full.niml.dset"

		# first, generate a group mask - the nodes that are common to all subjects

		# gather a list of the surfaces
		surfs = []

		for subj_id in subj_ids:

			subj_conf = ns_aperture.config.get_conf( subj_id )
			subj_paths = ns_aperture.paths.get_subj_paths( subj_conf )

			subj_surf = subj_paths.summ.mask.full( hemi_ext )

			surfs.append( subj_surf )

		mask_surf = paths.surf_mask.full( hemi_ext )

		fmri_tools.utils.group_surf_mask( surf_paths = surfs,
		                                  mask_path = mask_surf
		                                )

		# then, back project to a volume

		rep_subj_conf = ns_aperture.config.get_conf( conf.ana.rep_subj_id )
		rep_subj_paths = ns_aperture.paths.get_subj_paths( rep_subj_conf )

		rep_spec = rep_subj_paths.reg.std_spec.full( "_" + hemi + ".spec" )
		rep_ref = rep_subj_paths.reg.anat_reg.full( "+orig" )
		rep_vol_mask = paths.rep_vol_mask.full( "_" + hemi + "+orig" )

		cmd = [ "3dSurf2Vol",
		        "-spec", rep_spec,
		        "-surf_A", "smoothwm",
		        "-surf_B", "pial",
		        "-sdata", mask_surf,
		        "-grid_parent", rep_subj_paths.reg.mean.full( "+orig" ),
		        "-sv", rep_ref,
		        "-map_func", "max",
		        "-prefix", rep_vol_mask,
		        "-overwrite"
		      ]

		fmri_tools.utils.run_cmd( " ".join( cmd ) )


def coh_test( conf, paths ):
	"Tests the coherence GLM"

	os.chdir( paths.coh_test.dir() )

	subj_ids = conf.all_subj.subj.keys()

	os.chdir( paths.base.dir() )

	for hemi in [ "lh", "rh" ]:

		hemi_ext = "-std_{h:s}".format( h = hemi )

		# brick in the GLM file corresponding to the coherent regressor
		beta_brick = "3"

		cmd = [ "3dttest++", "-setA" ]

		subj_ids = conf.all_subj.subj.keys()

		for subj_id in subj_ids:

			subj_conf = ns_aperture.config.get_conf( subj_id )
			subj_paths = ns_aperture.paths.get_subj_paths( subj_conf )

			beta_path = subj_paths.ana.glm.full( hemi_ext + "-full.niml.dset" )

			cmd.append( "'" + beta_path + "[" + beta_brick + "]'" )

		cmd.extend( [ "-prefix", paths.coh_test.file( hemi_ext + "-full.niml.dset" ),
		              "-mask", paths.surf_mask.file( hemi_ext + "-full.niml.dset" ),
		              "-overwrite"
		            ]
		          )

		fmri_tools.utils.run_cmd( " ".join( cmd ) )


def cluster_sim( conf, paths ):
	"Generates a FWE cluster threshold"

	for hemi in [ "lh", "rh" ]:

		script_path = paths.clust_script.full( "_" + hemi + ".sh" )

		rep_subj_conf = ns_aperture.config.get_conf( conf.ana.rep_subj_id )
		rep_subj_paths = ns_aperture.paths.get_subj_paths( rep_subj_conf )

		rep_spec = rep_subj_paths.reg.std_spec.full( "_" + hemi + ".spec" )
		rep_ref = rep_subj_paths.reg.anat_reg.full( "+orig" )
		rep_vol_mask = paths.rep_vol_mask.full( "_" + hemi + "+orig" )

		fmri_tools.utils.cluster_sim( spec_path = rep_spec,
		                              mask_path = rep_vol_mask,
		                              ref_path = rep_ref,
		                              results_dir = paths.clust_sim.full( "_" + hemi ),
		                              script_path = script_path,
		                              pvals = [ conf.ana.p_height_thr ],
		                              blur = conf.ana.smooth_fwhm,
		                              sigma = conf.ana.smooth_sigma,
		                              stop_after_script = True,
		                              verbosity = 2
		                            )


def cluster( conf, paths, dt ):
	"Applies the cluster threshold"

	for hemi in [ "lh", "rh" ]:

		hemi_ext = "-std_{h:s}".format( h = hemi )

		sim_path = os.path.join( paths.clust_sim.full( "_" + hemi ),
		                         "z.max.area." + str( conf.ana.p_height_thr )
		                       )

		min_area = fmri_tools.utils.get_clust_thresh( sim_path )

		spec_path = ( paths.avg / "SUMA" ).full( "std.141.ns_aperture_avg_" + hemi + ".spec" )

		if dt == "coh":

			surf_path = paths.coh_test.full( hemi_ext + "-full.niml.dset" )
			surf_brick = 0

			clust_base = paths.coh_clust.full( hemi_ext + "-full.niml.dset" )
			thresh = 4.604
			abs_thresh = True
			thresh_brick = 1

			clust_base_path = paths.coh_clust

		elif dt == "loc":

			surf_path = paths.loc_test.full( hemi_ext + "-full.niml.dset" )
			surf_brick = 0

			clust_base = paths.loc_clust.full( hemi_ext + "-full.niml.dset" )
			thresh = 2.78
			abs_thresh = False
			thresh_brick = 1

			clust_base_path = paths.loc_clust

		elif dt == "acc":

			hemi_ext = "_" + hemi

			surf_path = paths.acc.full( "_" + hemi + ".niml.dset" )
			surf_brick = 0

			clust_base = paths.acc_clust.full( "_" + hemi + "-full.niml.dset" )
			thresh = 3.7469  # p < 0.02
			abs_thresh = False
			thresh_brick = 1

			clust_base_path = paths.acc_clust

		fmri_tools.utils.surf_cluster( surf_path = surf_path,
		                               surf_brick = surf_brick,
		                               spec_path = spec_path,
		                               clust_path = clust_base,
		                               thresh = thresh,
		                               abs_thresh = abs_thresh,
		                               thresh_brick = thresh_brick,
		                               min_area = min_area,
		                               ref_surf = "midway"
		                             )

		clust_ext = "-full_Clustered_e1_a" + str( min_area ) + ".niml.dset"
		clust_path = clust_base_path.full( hemi_ext + clust_ext )

		full_clust = clust_base_path.full( hemi_ext + "-full.niml.dset" )

		fmri_tools.utils.sparse_to_full( clust_path, full_clust, "ld141" )



def coh_effect_size( conf, paths ):
	"Calculates the effect size (avg beta / std beta ) of the coh cluster"

	subj_ids = conf.all_subj.subj.keys()

	data = np.empty( ( len( subj_ids ), 2 ) )
	data.fill( np.NAN )

	for ( i_subj, subj_id ) in enumerate( subj_ids ):

		subj_conf = ns_aperture.config.get_conf( subj_id )
		subj_paths = ns_aperture.paths.get_subj_paths( subj_conf )

		for ( i_hemi, hemi ) in enumerate(  [ "lh", "rh" ] ):

			beta_path = subj_paths.ana.clust.full( "-std_{h:s}.txt".format( h = hemi ) )

			subj_beta = np.loadtxt( beta_path )

			data[ i_subj, i_hemi ] = np.mean( subj_beta )

	data = np.mean( data, axis = 1 )

	mean = np.mean( data )
	sd = np.std( data, ddof = 1 )
	se = sd / np.sqrt( len( data ) )

	g1 = mean / sd


	print "Mean: " + str( mean )
	print "SE: " + str( se )
	print "g1: " + str( g1 )


def vis_loc_std( conf, paths ):
	"Analyses LOC and hMT localisers across subjects"

	subj_ids = conf.all_subj.subj.keys()

	os.chdir( paths.vis_loc.dir() )

	for hemi in [ "lh", "rh" ]:

		for dt in [ "obj", "motion" ]:

			glm_paths = []

			for subj_id in subj_ids:

				subj_conf = ns_aperture.config.get_conf( subj_id )
				subj_paths = ns_aperture.paths.get_subj_paths( subj_conf )

				# workaround for s1011, who had wedges collected in a separate session
				if ( subj_id == "s1011" ):
					vis_loc_date = "20111214"
				else:
					vis_loc_date = subj_conf.subj.vis_loc_date

				vis_loc_dir = os.path.join( "/labs/olmanlab/FsVisLoc",
				                            subj_id,
				                            vis_loc_date,
				                            "dt",
				                            dt
				                          )


				glm_path = os.path.join( vis_loc_dir,
				                         ( subj_id + "_vis_loc_" + vis_loc_date +
				                           "_" + dt + "-glm-std_" + hemi + ".niml.dset"
				                         )
				                       )

				glm_full_path = os.path.join( vis_loc_dir,
				                              ( subj_id + "_vis_loc_" + vis_loc_date +
				                                "_" + dt + "-glm-std_" + hemi + "-full.niml.dset"
				                              )
				                            )

				fmri_tools.utils.sparse_to_full( glm_path, glm_full_path, "ld141" )

				glm_paths.append( glm_full_path )


			cmd = [ "3dttest++", "-setA" ]

			cmd.extend( glm_paths )

			cmd.extend( [ "-prefix", paths.vis_loc.file( "-" + dt + "-std_" + hemi + "-full.niml.dset" ),
			              "-mask", paths.surf_mask.full( "-std_" + hemi + "-full.niml.dset" ),
			              "-overwrite"
			            ]
			          )

			fmri_tools.utils.run_cmd( " ".join( cmd ) )



def ret_std( conf, paths ):
	"Averages the ring and wedge maps across subjects"

	subj_ids = conf.all_subj.subj.keys()

	for dt in [ "wedge", "ring" ]:

		coef_paths = { "lh" : { "sin" : [], "cos" : [] },
		               "rh" : { "sin" : [], "cos" : [] }
		             }

		for subj_id in subj_ids:

			subj_conf = ns_aperture.config.get_conf( subj_id )
			subj_paths = ns_aperture.paths.get_subj_paths( subj_conf )

			# workaround for s1011, who had rings collected in a separate session
			if ( subj_id == "s1011" ) and ( dt == "ring" ):
				vis_loc_date = "20111214"
			else:
				vis_loc_date = subj_conf.subj.vis_loc_date


			vis_loc_dir = os.path.join( "/labs/olmanlab/FsVisLoc",
			                            subj_id,
			                            vis_loc_date,
			                            "dt",
			                            dt
			                          )

			for hemi in [ "lh", "rh" ]:


				ang_path = os.path.join( vis_loc_dir,
				                         ( subj_id +
				                           "_vis_loc_" +
				                           vis_loc_date + 
				                           "_" + dt + "-angle-std_" +
				                           hemi +
				                           ".niml.dset[0]"
				                         )
				                       )

				for coef in [ "sin", "cos" ]:

					out_ext = "-" + dt + "-" + coef + "-" + hemi + "-full.niml.dset"

					os.chdir( subj_paths.ana.coef.dir() )

					out_full_path = subj_paths.ana.coef.full( out_ext )
					out_file_path = subj_paths.ana.coef.file( out_ext )

					coef_paths[ hemi ][ coef ].append( out_full_path )

					cmd = [ "3dcalc",
					        "-a", ang_path,
					        "-expr", coef + "d(a)",
					        "-prefix", out_file_path,
					        "-pad_to_node", "ld141",
					        "-overwrite"
					      ]

					fmri_tools.utils.run_cmd( " ".join( cmd ) )

		os.chdir( paths.coef.dir() )

		# now to average across subjects
		letters = string.lowercase

		for hemi in [ "lh", "rh" ]:

			for coef in [ "sin", "cos" ]:

				cmd = [ "3dcalc" ]

				for ( letter, in_path ) in zip( letters, coef_paths[ hemi ][ coef ] ):
					cmd.extend( [ "-" + letter, in_path ] )

				cmd.extend( [ "-expr",
				              "mean(" + ",".join( letters[ :len( subj_ids ) ] ) + ")"
				            ]
				          )

				cmd.extend( [ "-prefix",
				              paths.coef.file( "-" + dt + "-" +
				                               coef +
				                               "-" +
				                               hemi +
				                               "-full.niml.dset"
				                             )
				            ]
				          )

				cmd.append( "-overwrite" )

				fmri_tools.utils.run_cmd( " ".join( cmd ) )

		for hemi in [ "lh", "rh" ]:

			# convert back to degrees
			cmd = [ "3dcalc",
			        "-a", paths.coef.file( "-" + dt + "-sin-" + hemi + "-full.niml.dset" ),
			        "-b", paths.coef.file( "-" + dt + "-cos-" + hemi + "-full.niml.dset" ),
			        "-expr", "'mod(mod(atan2(-a,b)*(180/PI)+360,360),360)'",
			        "-prefix", paths.coef.file( "-" + dt + "-" + hemi + "-full.niml.dset" ),
			        "-overwrite"
			      ]

			fmri_tools.utils.run_cmd( " ".join( cmd ) )


		for hemi in [ "lh", "rh" ]:

			# average the SNRs
			cmd = [ "3dcalc" ]

			for ( letter, subj_id ) in zip( letters, subj_ids ):

				subj_conf = ns_aperture.config.get_conf( subj_id )
				subj_paths = ns_aperture.paths.get_subj_paths( subj_conf )

				# workaround for s1011, who had rings collected in a separate session
				if ( subj_id == "s1011" ) and ( dt == "ring" ):
					vis_loc_date = "20111214"
				else:
					vis_loc_date = subj_conf.subj.vis_loc_date

				vis_loc_dir = os.path.join( "/labs/olmanlab/FsVisLoc",
				                            subj_id,
				                            vis_loc_date,
				                            "dt",
				                            dt
				                          )


				snr_path = os.path.join( vis_loc_dir,
				                         ( subj_id +
				                           "_vis_loc_" +
				                           vis_loc_date +
				                           "_" + dt + "-angle-std_" +
				                           hemi +
				                           ".niml.dset[1]"
				                         )
				                       )

				out_path = subj_paths.ana.snr.full( "-" + dt + "-" + hemi + "-full.niml.dset" )

				fmri_tools.utils.sparse_to_full( in_dset = snr_path,
				                                 out_dset = out_path,
				                                 pad_node = "ld141"
				                               )

				cmd.extend( [ "-" + letter, out_path ] )

			cmd.extend( [ "-expr",
			              "mean(" + ",".join( letters[ :len( subj_ids ) ] ) + ")"
			            ]
			          )

			cmd.extend( [ "-prefix",
			              paths.coef.file( "-" + dt + "-snr-" + hemi + "-full.niml.dset" )
			            ]
			          )

			cmd.append( "-overwrite" )

			fmri_tools.utils.run_cmd( " ".join( cmd ) )

			# combine the angle and snr
			cmd = [ "3dbucket",
			        "-prefix", paths.angle.file( "-" + dt + "-" + hemi + "-full.niml.dset" ),
			        "-overwrite",
			        paths.coef.file( "-" + dt + "-" + hemi + "-full.niml.dset" ),
			        paths.coef.file( "-" + dt + "-snr-" + hemi + "-full.niml.dset" )
			      ]

			fmri_tools.utils.run_cmd( " ".join( cmd ) )


def loc_test( conf, paths ):
	"Tests the localiser GLM"

	os.chdir( paths.loc.dir() )

	subj_ids = conf.all_subj.subj.keys()

	for hemi in [ "lh", "rh" ]:

		hemi_ext = "-std_{h:s}".format( h = hemi )

		# brick in the beta file corresponding to the all>0 regressor
		beta_brick = "5"

		cmd = [ "3dttest++", "-setA" ]

		subj_ids = conf.all_subj.subj.keys()

		for subj_id in subj_ids:

			subj_conf = ns_aperture.config.get_conf( subj_id )
			subj_paths = ns_aperture.paths.get_subj_paths( subj_conf )

			beta_path = subj_paths.loc.glm.full( hemi_ext + "-full.niml.dset" )

			cmd.append( "'" + beta_path + "[" + beta_brick + "]'" )

		cmd.extend( [ "-prefix", paths.loc_test.file( hemi_ext + "-full.niml.dset" ),
		              "-mask", paths.surf_mask.full( hemi_ext + "-full.niml.dset" ),
		              "-overwrite"
		            ]
		          )

		fmri_tools.utils.run_cmd( " ".join( cmd ) )


def mvpa_node_prep( conf, paths ):
	"Prepare the searchlight info for each node"

	rep_subj_conf = ns_aperture.config.get_conf( conf.ana.rep_subj_id )
	rep_subj_paths = ns_aperture.paths.get_subj_paths( rep_subj_conf )

	for hemi in [ "lh", "rh" ]:

		seed_nodes = np.loadtxt( rep_subj_paths.mvpa.nodes.full( "_" + hemi + ".txt" ) )

		surf_path = ( paths.avg / "SUMA" ).full( "std.141." + hemi + ".midway.asc" )

		pbar = progressbar.ProgressBar( widgets = [ progressbar.Percentage(),
		                                            progressbar.Bar()
		                                          ],
		                                maxval = seed_nodes.shape[ 0 ]
		                              ).start()

		with open( paths.sl_info.full( "_" + hemi + ".txt" ), "w" ) as node_file:

			for ( i_node, seed_node ) in enumerate( seed_nodes ):

				pbar.update( i_node )

				# save the seed node to a file
				np.savetxt( paths.sl_seed.full( ".txt" ), [ seed_node ], "%d" )

				cmd = [ "ROIgrow",
				        "-i", surf_path,
				        "-roi_nodes", paths.sl_seed.full( ".txt" ),
				        "-lim", "{r:.2f}".format( r = conf.ana.sl_r ),
				        "-prefix", paths.sl_disk.full( ".1D" ),
				        "-overwrite"
				      ]

				fmri_tools.utils.run_cmd( " ".join( cmd ), log_stdout = False )

				sl_nodes = np.loadtxt( paths.sl_disk.full( ".1D" ) )

				sl_nodes = np.concatenate( ( [ seed_node ], sl_nodes ) )
				sl_nodes = np.unique( sl_nodes )

				seed_str = ( "\t".join( [ "{n:.0f}".format( n = sl_node )
				                          for sl_node in sl_nodes
				                        ]
				                      ) +
				             "\n"
				           )

				node_file.write( seed_str )

		pbar.finish()


def mvpa_test( conf, paths ):
	"Tests the mvpa accuracies"

	os.chdir( paths.coh_test.dir() )

	subj_ids = conf.all_subj.subj.keys()

	os.chdir( paths.base.dir() )

	for hemi in [ "lh", "rh" ]:

		continue
		hemi_ext = "-std_{h:s}".format( h = hemi )

		test_paths = []

		subj_ids = conf.all_subj.subj.keys()

		for subj_id in subj_ids:

			subj_conf = ns_aperture.config.get_conf( subj_id )
			subj_paths = ns_aperture.paths.get_subj_paths( subj_conf )

			orig_mvpa_path = subj_paths.mvpa.acc.full( "_" + hemi + ".niml.dset" )
			mvpa_path = subj_paths.mvpa.acc.full( "-smooth_" + hemi + ".niml.dset" )
			spec_file = subj_paths.reg.std_spec.full( "_" + hemi + ".spec" )

			blur_path = subj_paths.func.surfs[ 0 ].full( "-std_" + hemi + ".niml.dset" )

			fmri_tools.preproc.surf_smooth( in_surf = orig_mvpa_path,
			                                out_surf = mvpa_path,
			                                spec_path = spec_file,
			                                target_fwhm = conf.ana.smooth_fwhm,
			                                surf = "std.141." + hemi + ".midway.asc",
			                                blur_master = blur_path,
			                                extra_params = [ "-sigma",
			                                                 conf.ana.smooth_sigma
			                                               ]
			                              )

			# now to subtract 50
			os.chdir( subj_paths.mvpa.acc.dir() )

			rel_acc_path = subj_paths.mvpa.acc.full( "-smooth_rel_" + hemi + ".niml.dset" )

			cmd = [ "3dcalc",
			        "-a", mvpa_path,
			        "-expr", "a-50",
			        "-prefix", rel_acc_path,
			        "-overwrite"
			      ]

			fmri_tools.utils.run_cmd( " ".join( cmd ) )

			test_paths.append( rel_acc_path )

		# t-test

		os.chdir( paths.acc.dir() )

		cmd = [ "3dttest++", "-setA" ]

		cmd.extend( test_paths )

		cmd.extend( [ "-prefix", paths.acc.full( "_" + hemi + ".niml.dset" ) ] )

		cmd.append( "-overwrite" )

		fmri_tools.utils.run_cmd( " ".join( cmd ) )

	cluster( conf, paths, "acc" )
