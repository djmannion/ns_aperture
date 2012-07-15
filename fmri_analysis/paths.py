"""Paths for the natural scenes apertures fMRI experiment"""

from __future__ import division

import os

import fmri_tools.paths

def _get_study_paths():
	"""Get the path structure for the study"""

	base_dir = "/labs/olmanlab/Data7T/NatSceneAperture/"

	study_paths = fmri_tools.paths.get_study_paths( base_dir )

	return study_paths


def _get_func_paths( conf, paths ):
	"Returns the paths associated with the functional acquisitions"

	func_dir = os.path.join( paths[ "study" ][ "subj_dir" ],
	                         conf[ "subj" ][ "subj_id" ],
	                         "func"
	                       )

	func_exp_dir = os.path.join( func_dir, "exp" )

	func_paths = fmri_tools.paths.get_func_paths( func_exp_dir,
	                                              conf[ "subj" ][ "subj_id" ],
	                                              conf[ "subj" ][ "n_runs" ],
	                                              conf[ "exp" ][ "id" ]
	                                            )

	func_paths[ "surf_files" ] = [ orig_file.replace( "orig", "surf" )
	                               for orig_file in func_paths[ "orig_files" ]
	                             ]

	paths[ "func_exp" ] = func_paths

	func_loc_dir = os.path.join( func_dir, "loc" )

	loc_paths = fmri_tools.paths.get_func_paths( func_loc_dir,
	                                             conf[ "subj" ][ "subj_id" ],
	                                             conf[ "subj" ][ "n_loc_runs" ],
	                                             "%s_loc" % conf[ "exp" ][ "id" ]
	                                           )

	paths[ "func_loc" ] = loc_paths

	return paths


def _get_summ_paths( conf, paths ):
	"Returns the paths for images summarising the functional session"

	summ_dir = os.path.join( paths[ "study" ][ "subj_dir" ],
	                         conf[ "subj" ][ "subj_id" ],
	                         "func"
	                       )

	summ_paths = fmri_tools.paths.get_func_summ_paths( summ_dir,
	                                                   conf[ "subj" ][ "subj_id" ],
	                                                   conf[ "exp" ][ "id" ]
	                                                 )

	log_file = "%s_%s-log.log" % ( conf[ "subj" ][ "subj_id" ],
	                               conf[ "exp" ][ "id" ]
	                             )

	summ_paths[ "log_file" ] = os.path.join( paths[ "study" ][ "subj_dir" ],
	                                         log_file
	                                       )

	paths[ "summ" ] = summ_paths

	return paths


def _get_fmap_paths( conf, paths ):
	"Returns the paths for the session's fieldmaps"

	fmap_dir = os.path.join( paths[ "study" ][ "subj_dir" ],
	                         conf[ "subj" ][ "subj_id" ],
	                         "fmap"
	                       )

	fmap_paths = fmri_tools.paths.get_fmap_paths( fmap_dir,
	                                              conf[ "subj" ][ "subj_id" ],
	                                              conf[ "exp" ][ "id" ],
	                                              conf[ "subj" ][ "n_fmaps" ]
	                                            )

	paths[ "fmap" ] = fmap_paths

	return paths


def _get_reg_paths( conf, paths ):
	"Returns the paths for the session's registration with FreeSurfer anatomy"

	reg_dir = os.path.join( paths[ "study" ][ "subj_dir" ],
	                        conf[ "subj" ][ "subj_id" ],
	                        "reg"
	                      )

	reg_paths = {}

	subj_id = conf[ "subj" ][ "subj_id" ]

	reg_paths[ "exp_anat" ] = os.path.join( reg_dir,
	                                        "%s_anat.nii" % subj_id
	                                      )

	reg_paths[ "rs_exp_anat" ] = os.path.join( reg_dir,
	                                           "%s_anat_rs+orig" % subj_id
	                                         )

	surf_dir = os.path.join( os.environ[ "SUBJECTS_DIR" ],
	                         subj_id,
	                         "SUMA"
	                       )

	reg_paths[ "surf_anat" ] = os.path.join( surf_dir,
	                                         "%s_SurfVol+orig" % subj_id
	                                       )

	reg_paths[ "reg" ] = os.path.join( reg_dir,
	                                   "%s_SurfVol_Alnd_Exp+orig" % subj_id
	                                 )

	reg_paths[ "spec" ] = os.path.join( surf_dir,
	                                    "%s_" % subj_id
	                                  )

	paths[ "reg" ] = reg_paths

	return paths


def _get_log_paths( conf, paths ):
	"""Get the paths for the logfiles"""

	subj_id = conf[ "subj" ][ "subj_id" ]

	log = {}

	log_dir = os.path.join( paths[ "study" ][ "subj_dir" ],
	                        subj_id,
	                        "logs"
	                      )

	log[ "base_dir" ] = log_dir

	log[ "design" ] = os.path.join( log_dir, "design.npy" )

	log[ "loc_design" ] = os.path.join( log_dir, "loc_design.npy" )

	log[ "seq_base" ] = os.path.join( log_dir,
	                                  "%s_ns_aperture_fmri_seq_" % subj_id
	                                )

	log[ "task_base" ] = os.path.join( log_dir,
	                                   "%s_ns_aperture_fmri_task_" % subj_id
	                                 )

	paths[ "log" ] = log

	return paths


def get_subj_paths( conf ):
	"""Get the path structure for a given subject"""

	paths = {}

	paths[ "study" ] = _get_study_paths()

	paths = _get_func_paths( conf, paths )

	paths = _get_summ_paths( conf, paths )

	paths = _get_fmap_paths( conf, paths )

	paths = _get_reg_paths( conf, paths )

	paths = _get_log_paths( conf, paths )

	return paths
