from __future__ import division

import os, os.path
import logging

import fmri_tools.preproc, fmri_tools.utils


def convert( conf, paths ):
	"""Converts the functionals and fieldmaps from dicom to nifti"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running conversion..." )

	# aggregate the dicom directories
	raw_paths = ( paths.func.raws +
	              [ paths.fmap.mag_raw, paths.fmap.ph_raw ]
	            )

	# aggregate the image paths
	img_paths = ( paths.func.origs +
	              [ paths.fmap.mag, paths.fmap.ph ]
	            )

	for ( raw_path, img_path ) in zip( raw_paths, img_paths ):

		fmri_tools.preproc.dcm_to_nii( dcm_path = raw_path.full(),
		                               nii_path = img_path.full( ".nii" )
		                             )

		# reorient
		fmri_tools.preproc.reorient_img( in_path = img_path.full( ".nii" ),
		                                 out_path = img_path.full( ".nii" ),
		                                 reorient_str = conf.acq.reshape_to_RAS
		                               )

	# check that they are all unique
	are_unique = fmri_tools.utils.files_are_unique( [ img_path.full( ".nii" )
	                                                  for img_path in img_paths
	                                                ]
	                                              )

	assert are_unique

	# files to go into the summary
	summ_paths = [ orig.full() for orig in paths.func.origs ]

	# make a summary image from the files
	fmri_tools.preproc.gen_sess_summ_img( epi_paths = summ_paths,
	                                      out_path = paths.summ.orig.full(),
	                                    )


def fieldmaps( conf, paths ):
	"""Prepare the fieldmaps"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running fieldmap preparation..." )

	epi_run = conf.subj.mot_base - 1

	# set a corrected EPI to define the space to resample to
	ref_epi = paths.func.origs[ epi_run ].full( ".nii[0]" )

	# want to calculate a coarse brain mask from the epi
	mask_path = paths.fmap.mask.full( ".nii" )
	epi_path = paths.summ.orig.full( ".nii[{n:d}]".format( n = epi_run ) )

	mask_cmd = [ "3dAutomask",
	             "-SI", "{n:d}".format( n = conf.subj.mask_SI ),
	             "-overwrite",
	             "-prefix", mask_path,
	             epi_path
	           ]

	fmri_tools.utils.run_cmd( " ".join( mask_cmd ) )

	fmri_tools.preproc.make_fieldmap( mag_path = paths.fmap.mag.full(),
	                                  ph_path = paths.fmap.ph.full(),
	                                  fmap_path = paths.fmap.fmap.full(),
	                                  delta_te_ms = conf.acq.delta_te_ms,
	                                  ref_img = ref_epi,
	                                  recentre_ph = "mean",
	                                  recentre_mask = mask_path,
	                                  strip_params = [ "-surface_coil" ],
	                                  strip_mag = False
	                                )


def mc_unwarp( conf, paths ):
	"""Combined motion and distortion correction"""

	logger = logging.getLogger( __name__ )
	logger.info( "Running motion and distortion correction..." )

	orig_paths = [ orig_path.full() for orig_path in paths.func.origs ]
	corr_paths = [ corr_path.full() for corr_path in paths.func.corrs ]
	uw_paths = [ uw_path.full() for uw_path in paths.func.uws ]

	fmap_path = paths.fmap.fmap.full()

	fmri_tools.preproc.mc_unwarp( orig_paths = orig_paths,
	                              corr_paths = corr_paths,
	                              uw_paths = uw_paths,
	                              fmap_path = fmap_path,
	                              dwell_ms = conf.acq.dwell_ms,
	                              uw_direction_fsl = conf.acq.ph_encode_dir,
	                              uw_direction_afni = conf.acq.ph_enc_afni,
	                              voxel_size = conf.acq.vox_size,
	                              i_vol_base = conf.acq.vol_base,
	                              i_run_base = conf.subj.mot_base - 1,
	                              fugue_params = conf.acq.fugue_params
	                            )

	fmri_tools.preproc.gen_sess_summ_img( epi_paths = uw_paths,
	                                      out_path = paths.summ.uw.full()
	                                    )

	fmri_tools.preproc.mean_image( uw_paths, paths.summ.mean.full() )


def sess_reg( conf, paths ):
	"""Coregisters the session anatomical with its mean functional"""

	logger = logging.getLogger( __name__ )

	logger.info( "Running registration..." )

	if conf.subj.extra_al_params is not None:
		extra_al_params = conf.subj.extra_al_params
	else:
		extra_al_params = None

	fmri_tools.preproc.img_reg( reg_dir = paths.reg.base.full(),
	                            base_file = paths.reg.mean.file( "+orig" ),
	                            mov_file = paths.reg.anat_ref.file( "+orig" ),
	                            extra_al_params = extra_al_params,
	                            epi_strip = "None" #"3dAutomask"
	                          )



def vol_to_surf( conf, paths, std_surf = False ):
	"""Converts the functional volume-based images to SUMA surfaces."""

	logger = logging.getLogger( __name__ )
	logger.info( "Running volume to surface projection..." )

	start_dir = os.getcwd()

	for ( uw_file, surf_file, run_dir ) in zip( paths.func.uws,
	                                            paths.func.surfs,
	                                            paths.func.runs
	                                          ):

		os.chdir( run_dir.full() )

		for hemi in [ "lh", "rh" ]:

			if std_surf:
				spec_file = paths.reg.std_spec.full( "_{hemi:s}.spec".format( hemi = hemi ) )
				surf_path = surf_file.full( "-std_{h:s}.niml.dset".format( h = hemi ) )
			else:
				spec_file = paths.reg.spec.full( "_{hemi:s}.spec".format( hemi = hemi ) )
				surf_path = surf_file.full( "_{h:s}.niml.dset".format( h = hemi ) )

			surf_cmd = [ "3dVol2Surf",
			             "-spec", spec_file,
			             "-surf_A", "smoothwm",
			             "-surf_B", "pial",
			             "-map_func", "ave",
			             "-f_steps", "15",
			             "-f_index", "nodes",
			             "-sv", paths.reg.anat_reg.full( "+orig" ),
			             "-grid_parent", uw_file.full( ".nii" ),
			             "-out_niml", surf_path,
			             "-overwrite"
			           ]

			fmri_tools.utils.run_cmd( " ".join( surf_cmd ) )

	os.chdir( start_dir )


def smooth_surfs( conf, paths ):
	"Smooth the surfaces for group analysis"

	logger = logging.getLogger( __name__ )
	logger.info( "Smoothing group surfaces..." )

	start_dir = os.getcwd()

	for run_num in xrange( conf.subj.n_runs ):

		surf_file = paths.func.surfs[ run_num - 1 ]
		run_dir = paths.func.runs[ run_num - 1 ]

		os.chdir( run_dir.full() )

		for hemi in [ "lh", "rh" ]:

			spec_file = paths.reg.std_spec.full( "_{hemi:s}.spec".format( hemi = hemi ) )
			surf_path = surf_file.full( "-std_{h:s}.niml.dset".format( h = hemi ) )
			smooth_path = surf_file.full( "-smooth-std_{h:s}.niml.dset".format( h = hemi ) )

			fmri_tools.preproc.surf_smooth( in_surf = surf_path,
			                                out_surf = smooth_path,
			                                spec_path = spec_file,
			                                target_fwhm = conf.ana.smooth_fwhm,
			                                surf_A = "std.141." + hemi + ".smoothwm.asc",
			                                surf_B = "std.141." + hemi + ".pial.asc",
			                                extra_params = [ "-sigma",
			                                                 conf.ana.smooth_sigma
			                                               ]
			                              )

	os.chdir( start_dir )


def filter_tc( conf, paths ):
	"Filter the surface timeseries for MVPA analyses"

	logger = logging.getLogger( __name__ )
	logger.info( "Filtering timeseries..." )

	# minus 1 because exp_runs is 1 based
	surfs = [ paths.func.surfs[ i_surf - 1 ] for i_surf in conf.subj.exp_runs ]
	filts = [ paths.func.filts[ i_surf - 1 ] for i_surf in conf.subj.exp_runs ]

	for hemi in [ "lh", "rh" ]:

		for ( surf, filt ) in zip( surfs, filts ):

			os.chdir( filt.dir() )

			surf_path = surf.file( "-std_{h:s}.niml.dset".format( h = hemi ) )
			filt_path = filt.file( "-std_{h:s}.niml.dset".format( h = hemi ) )

			filt_cmd = [ "3dDetrend",
			             "-prefix", filt_path,
			             "-polort", conf.ana.mvpa_filt_ord,
			             "-overwrite",
			             surf_path
			           ]

			fmri_tools.utils.run_cmd( " ".join( filt_cmd ) )

			# pad to full because it needs to be combined with the group mask

			full_filt_path = filt.file( "-std_{h:s}-full.niml.dset".format( h = hemi ) )

			fmri_tools.utils.sparse_to_full( filt_path,
			                                 full_filt_path,
			                                 pad_node = "ld141"
			                               )


def surf_mask( conf, paths ):
	"Creates a mask of all nodes nonzero in all runs"

	logger = logging.getLogger( __name__ )
	logger.info( "Creating surface masks..." )

	for hemi in [ "lh", "rh" ]:

		surfs = [ surf.full( "-std_{h:s}.niml.dset[0]".format( h = hemi ) )
		          for surf in paths.func.surfs
		        ]

		mask_path = paths.summ.mask.full( "-std_{h:s}.niml.dset".format( h = hemi ) )

		fmri_tools.utils.group_surf_mask( surf_paths = surfs,
		                                  mask_path = mask_path
		                                )

		full_mask_path = paths.summ.mask.full( "-std_{h:s}-full.niml.dset".format( h = hemi ) )

		fmri_tools.utils.sparse_to_full( mask_path,
		                                 full_mask_path,
		                                 pad_node = "ld141"
		                               )

