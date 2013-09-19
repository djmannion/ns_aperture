
import os.path
import string

import fmri_tools.utils

import ns_aperture.config, ns_aperture.paths


def avg_phase_surfs( conf, paths ):

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
		        "-prefix", paths.wedge_angle.file( "-angle-" + hemi + "-full.niml.dset" ),
		        "-overwrite",
		        paths.wedge_coef.file( "-" + hemi + "-full.niml.dset" ),
		        paths.wedge_coef.file( "-snr-" + hemi + "-full.niml.dset" )
		      ]

		fmri_tools.utils.run_cmd( " ".join( cmd ) )