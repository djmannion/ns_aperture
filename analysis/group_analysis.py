
import os.path
import string

import fmri_tools.utils

import ns_aperture.config, ns_aperture.paths


def coh_test( conf, paths ):
	"Tests the coherence GLM"

	os.chdir( paths.coh_test.dir() )

	subj_ids = conf.all_subj.subj.keys()

	os.chdir( paths.base.dir() )

	for hemi in [ "lh", "rh" ]:

		hemi_ext = "-std_{h:s}".format( h = hemi )

		# first, make a mask from the overlapping regions
		cmd = [ "3dcalc" ]

		for ( letter, subj_id ) in zip( string.letters, subj_ids ):

			subj_conf = ns_aperture.config.get_conf( subj_id )
			subj_paths = ns_aperture.paths.get_subj_paths( subj_conf )

			cmd.extend( [ "-" + letter,
			              subj_paths.ana.glm.full( hemi_ext + "-full.niml.dset" )
			            ]
			          )

		expr = "1*and("
		expr += ",".join( [ "notzero(" + x + ")"
		                    for x in string.letters[ :len( subj_ids ) ]
		                  ]
		                )
		expr += ")"

		cmd.extend( [ "-expr", "'" + expr + "'" ] )

		cmd.extend( [ "-prefix",
		              paths.mask.file( hemi_ext + "-full.niml.dset" )
		            ]
		          )

		cmd.append( "-overwrite" )

		fmri_tools.utils.run_cmd( " ".join( cmd ) )

		# brick in the beta file corresponding to the coherent regressor
		beta_brick = "30"


		cmd = [ "3dttest++", "-setA" ]

		subj_ids = conf.all_subj.subj.keys()

		for subj_id in subj_ids:

			subj_conf = ns_aperture.config.get_conf( subj_id )
			subj_paths = ns_aperture.paths.get_subj_paths( subj_conf )

			beta_path = subj_paths.ana.beta.full( hemi_ext + "-full.niml.dset" )

			cmd.append( "'" + beta_path + "[" + beta_brick + "]'" )

		cmd.extend( [ "-prefix", paths.coh_test.file( hemi_ext + "-full.niml.dset" ),
		              "-mask", paths.mask.file( hemi_ext + "-full.niml.dset" ),
		              "-overwrite"
		            ]
		          )

		fmri_tools.utils.run_cmd( " ".join( cmd ) )


def cluster_sim( conf, paths ):
	"Generates a FWE cluster threshold"

	subj_ids = conf.all_subj.subj.keys()

	for hemi in [ "lh", "rh" ]:

		hemi_ext = "-std_" + hemi + "-full.niml.dset"

		# first, generate a group mask - the nodes that are common to all subjects

		# gather a list of the surfaces
		surfs = []

		for subj_id in subj_ids:

			subj_conf = ns_aperture.config.get_conf( subj_id )
			subj_paths = ns_aperture.paths.get_subj_paths( subj_conf )

			subj_surf = subj_paths.ana.glm.full( hemi_ext + "[0]" )

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

		# finally, run the cluster sim

		script_path = paths.clust_script.full( "_" + hemi + ".sh" )

		fmri_tools.utils.cluster_sim( spec_path = rep_spec,
		                              mask_path = rep_vol_mask,
		                              ref_path = rep_ref,
		                              results_dir = paths.clust_sim.full( "_" + hemi ),
		                              script_path = script_path,
		                              pvals = [ conf.ana.p_height_thr ],
		                              blur = conf.ana.smooth_fwhm,
		                              sigma = conf.ana.smooth_sigma,
		                              verbosity = 2
		                            )


def avg_phase_surfs( conf, paths ):
	"Averages the wedge maps across subjects"

	subj_ids = conf.all_subj.subj.keys()

	coef_paths = { "lh" : { "sin" : [], "cos" : [] },
	               "rh" : { "sin" : [], "cos" : [] }
	             }

	for subj_id in subj_ids:

		subj_conf = ns_aperture.config.get_conf( subj_id )
		subj_paths = ns_aperture.paths.get_subj_paths( subj_conf )

		vis_loc_dir = os.path.join( "/labs/olmanlab/FsVisLoc",
		                            subj_id,
		                            subj_conf.subj.vis_loc_date,
		                            "dt",
		                            "wedge"
		                          )

		for hemi in [ "lh", "rh" ]:

			ang_path = os.path.join( vis_loc_dir,
			                         ( subj_id +
			                           "_vis_loc_" +
			                           subj_conf.subj.vis_loc_date +
			                           "_wedge-angle-std_" +
			                           hemi +
			                           ".niml.dset[0]"
			                         )
			                       )

			for coef in [ "sin", "cos" ]:

				out_ext = "-" + coef + "-" + hemi + "-full.niml.dset"

				os.chdir( subj_paths.ana.wedge_coef.dir() )

				out_full_path = subj_paths.ana.wedge_coef.full( out_ext )
				out_file_path = subj_paths.ana.wedge_coef.file( out_ext )

				coef_paths[ hemi ][ coef ].append( out_full_path )

				cmd = [ "3dcalc",
				        "-a", ang_path,
				        "-expr", coef + "d(a)",
				        "-prefix", out_file_path,
				        "-pad_to_node", "ld141",
				        "-overwrite"
				      ]

				fmri_tools.utils.run_cmd( " ".join( cmd ) )

	os.chdir( paths.wedge_coef.dir() )

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
			              paths.wedge_coef.file( "-" +
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
		        "-a", paths.wedge_coef.file( "-sin-" + hemi + "-full.niml.dset" ),
		        "-b", paths.wedge_coef.file( "-cos-" + hemi + "-full.niml.dset" ),
		        "-expr", "'mod(mod(atan2(-a,b)*(180/PI)+360,360),360)'",
		        "-prefix", paths.wedge_coef.file( "-" + hemi + "-full.niml.dset" ),
		        "-overwrite"
		      ]

		fmri_tools.utils.run_cmd( " ".join( cmd ) )


	for hemi in [ "lh", "rh" ]:

		# average the SNRs
		cmd = [ "3dcalc" ]

		for ( letter, subj_id ) in zip( letters, subj_ids ):

			subj_conf = ns_aperture.config.get_conf( subj_id )
			subj_paths = ns_aperture.paths.get_subj_paths( subj_conf )

			vis_loc_dir = os.path.join( "/labs/olmanlab/FsVisLoc",
			                            subj_id,
			                            subj_conf.subj.vis_loc_date,
			                            "dt",
			                            "wedge"
			                          )


			snr_path = os.path.join( vis_loc_dir,
			                         ( subj_id +
			                           "_vis_loc_" +
			                           subj_conf.subj.vis_loc_date +
			                           "_wedge-angle-std_" +
			                           hemi +
			                           ".niml.dset[1]"
			                         )
			                       )

			out_path = subj_paths.ana.wedge_snr.full( "-" + hemi + "-full.niml.dset" )

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
		              paths.wedge_coef.file( "-snr-" + hemi + "-full.niml.dset" )
		            ]
		          )

		cmd.append( "-overwrite" )

		fmri_tools.utils.run_cmd( " ".join( cmd ) )

		# combine the angle and snr
		cmd = [ "3dbucket",
		        "-prefix", paths.wedge_angle.file( "-" + hemi + "-full.niml.dset" ),
		        "-overwrite",
		        paths.wedge_coef.file( "-" + hemi + "-full.niml.dset" ),
		        paths.wedge_coef.file( "-snr-" + hemi + "-full.niml.dset" )
		      ]

		fmri_tools.utils.run_cmd( " ".join( cmd ) )
