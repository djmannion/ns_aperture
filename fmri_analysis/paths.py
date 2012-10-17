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

	func_paths = fmri_tools.paths.get_func_paths( func_dir,
	                                              conf[ "subj" ][ "subj_id" ],
	                                              conf[ "exp" ][ "id" ],
	                                              conf[ "subj" ][ "n_runs" ],
	                                            )

	paths[ "func" ] = func_paths

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

	summ_paths[ "mot_est_file" ] = summ_paths[ "mot_est_file" ].replace( "npy",
	                                                                     "txt"
	                                                                   )

	log_file = "%s_%s-log.log" % ( conf[ "subj" ][ "subj_id" ],
	                               conf[ "exp" ][ "id" ]
	                             )

	summ_paths[ "log_file" ] = os.path.join( paths[ "study" ][ "subj_dir" ],
	                                         conf[ "subj" ][ "subj_id" ],
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

	reg_paths = fmri_tools.paths.get_reg_paths( reg_dir,
	                                            conf[ "subj" ][ "subj_id" ],
	                                            conf[ "exp" ][ "id" ]
	                                          )

	reg_paths[ "reg_file" ] = reg_paths[ "algn_file" ]

	paths[ "reg" ] = reg_paths

	return paths


def _get_task_paths( conf, paths ):
	"""Get the paths for the task"""

	subj_id = conf[ "subj" ][ "subj_id" ]

	task = {}

	task_dir = os.path.join( paths[ "study" ][ "subj_dir" ],
	                        subj_id,
	                        "task"
	                      )

	task[ "base_dir" ] = task_dir

	task[ "resp" ] = os.path.join( task_dir,
	                               "%s_ns_aperture_fmri_resp.npy" % subj_id
	                             )

	task[ "hist" ] = os.path.join( task_dir,
	                               "%s_ns_aperture_fmri_resp_hist.npy" % subj_id
	                             )

	paths[ "task" ] = task

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

	log[ "seq_base" ] = os.path.join( log_dir,
	                                  "%s_ns_aperture_fmri_seq_" % subj_id
	                                )

	log[ "task_base" ] = os.path.join( log_dir,
	                                   "%s_ns_aperture_fmri_task_" % subj_id
	                                 )

	paths[ "log" ] = log

	return paths


def _get_ana_paths( conf, paths ):
	"""Get the paths for the analysis"""

	subj_id = conf[ "subj" ][ "subj_id" ]
	exp_id = conf[ "exp" ][ "id" ]

	ana = {}

	ana_dir = os.path.join( paths[ "study" ][ "subj_dir" ],
	                        subj_id,
	                        "analysis"
	                      )

	ana[ "base_dir" ] = ana_dir

	ana[ "time_files" ] = os.path.join( ana_dir,
	                                    "%s_%s-stim_times_" % ( subj_id, exp_id )
	                                  )

	ana[ "bl_poly" ] = os.path.join( ana[ "base_dir" ],
	                                 "%s_%s-bl_poly.txt" % ( subj_id, exp_id )
	                               )

	ana[ "cens" ] = os.path.join( ana[ "base_dir" ],
	                              "%s_%s-cens.txt" % ( subj_id, exp_id )
	                            )

	ana[ "glm" ] = os.path.join( ana[ "base_dir" ],
	                             "%s_%s-glm" % ( subj_id, exp_id )
	                           )

	ana[ "beta" ] = os.path.join( ana[ "base_dir" ],
	                              "%s_%s-beta" % ( subj_id, exp_id )
	                            )

	ana[ "bltc" ] = os.path.join( ana[ "base_dir" ],
	                              "%s_%s-bltc" % ( subj_id, exp_id )
	                            )

	ana[ "bl" ] = os.path.join( ana[ "base_dir" ],
	                            "%s_%s-bl" % ( subj_id, exp_id )
	                          )

	ana[ "psc" ] = os.path.join( ana[ "base_dir" ],
	                             "%s_%s-psc" % ( subj_id, exp_id )
	                           )

	paths[ "ana" ] = ana

	return paths


def _get_loc_paths( conf, paths ):
	"""Get the paths for the localiser"""

	subj_id = conf[ "subj" ][ "subj_id" ]
	id = conf[ "exp" ][ "id_loc" ]

	loc = {}

	loc_dir = os.path.join( paths[ "study" ][ "subj_dir" ],
	                        subj_id,
	                        "loc"
	                      )

	loc[ "base_dir" ] = loc_dir

	loc[ "time_files" ] = os.path.join( loc_dir,
	                                    "%s_%s-stim_times_" % ( subj_id, id )
	                                  )

	loc[ "glm" ] = os.path.join( loc_dir,
	                             "%s_%s-glm" % ( subj_id, id )
	                           )

	loc[ "roi_stat" ] = os.path.join( loc_dir,
	                                  "%s_%s-roi_stat" % ( subj_id, id )
	                                )

	loc[ "parc" ] = os.path.join( loc_dir,
	                              "%s_%s-parc" % ( subj_id, id )
	                            )

	loc[ "roi_parc" ] = os.path.join( loc_dir,
	                                  "%s_%s-roi_parc" % ( subj_id, id )
	                                )

	paths[ "loc" ] = loc

	return paths


def _get_roi_paths( conf, paths ):
	"""Get the paths for the ROI data"""

	subj_id = conf[ "subj" ][ "subj_id" ]
	exp_id = conf[ "exp" ][ "id" ]

	rois = {}

	rois[ "base_dir" ] = os.path.join( paths[ "study" ][ "subj_dir" ],
	                                   subj_id,
	                                   "rois"
	                                 )

	rois[ "dset" ] = os.path.join( rois[ "base_dir" ],
	                               "%s_%s-rois" % ( subj_id, exp_id )
	                             )

	rois[ "fov_dset" ] = os.path.join( rois[ "base_dir" ],
	                                   "%s_%s-fov" % ( subj_id, exp_id )
	                                 )

	rois[ "psc" ] = os.path.join( rois[ "base_dir" ],
	                              "%s_%s-psc" % ( subj_id, exp_id )
	                            )

	rois[ "parc_psc" ] = os.path.join( rois[ "base_dir" ],
	                                   "%s_%s-parc_psc" % ( subj_id, exp_id )
	                                 )

	rois[ "raw_adj_tc" ] = os.path.join( rois[ "base_dir" ],
	                                     "%s_%s-raw_adj_tc" % ( subj_id, exp_id )
	                                   )

	rois[ "pred_adj_tc" ] = os.path.join( rois[ "base_dir" ],
	                                      "%s_%s-pred_adj_tc" % ( subj_id, exp_id )
	                                    )

	paths[ "rois" ] = rois

	return paths


def _get_mvpa_paths( conf, paths ):
	"""Get the paths for the MVPA analysis"""

	subj_id = conf[ "subj" ][ "subj_id" ]
	id = "ns_aperture_mvpa"

	mvpa = {}

	mvpa_dir = os.path.join( paths[ "study" ][ "subj_dir" ],
	                         subj_id,
	                         "mvpa"
	                       )

	mvpa[ "base_dir" ] = mvpa_dir


	mvpa[ "rfe_base_dir" ] = os.path.join( mvpa_dir, "rfe" )

	filt_files = [ os.path.join( mvpa_dir,
	                             "filt",
	                             os.path.split( surf_file )[ 1 ].replace( "surf", "surf_filt" )
	                           )
	               for surf_file in paths[ "func" ][ "surf_files" ]
	             ]

	mvpa[ "filt_files" ] = filt_files[ :conf[ "subj" ][ "n_exp_runs" ] ]

	mvpa[ "data" ] = os.path.join( mvpa_dir,
	                               "data",
	                               "%s_%s-data" % ( subj_id, id )
	                             )

	mvpa[ "parc_data" ] = os.path.join( mvpa_dir,
	                                    "data",
	                                    "%s_%s-parc_data" % ( subj_id, id )
	                                  )

	mvpa[ "parc_info" ] = os.path.join( mvpa_dir,
	                                    "data",
	                                    "%s_%s-parc_info" % ( subj_id, id )
	                                  )

	mvpa[ "data_info" ] = os.path.join( mvpa_dir,
	                                    "data",
	                                    "%s_%s-data_info" % ( subj_id, id )
	                                  )

	mvpa[ "blk_data_info" ] = os.path.join( mvpa_dir,
	                                        "data",
	                                        "%s_%s-blk_data_info.npy" % ( subj_id, id )
	                                      )

	mvpa[ "blk_data" ] = os.path.join( mvpa_dir,
	                                   "data",
	                                   "%s_%s-blk_parc_data" % ( subj_id, id )
	                                 )

	mvpa[ "nodes" ] = os.path.join( mvpa_dir,
	                                "%s_%s-nodes" % ( subj_id, id )
	                              )

	mvpa[ "acc_base" ] = os.path.join( "%s_%s-acc" % ( subj_id, id ) )

	mvpa[ "node_wgt_base" ] = os.path.join( "%s_%s-node_wgt" % ( subj_id, id ) )

	mvpa[ "wgt_base" ] = os.path.join( "%s_%s-wgt" % ( subj_id, id ) )

	mvpa[ "acc_mean" ] = os.path.join( mvpa_dir,
	                                   "%s_%s-acc_mean" % ( subj_id, id )
	                                 )

	mvpa[ "log_file" ] = os.path.join( mvpa_dir,
	                                   "%s_%s-mvpa_log.txt" % ( subj_id, id )
	                                 )

	mvpa[ "wgt_fov" ] = os.path.join( mvpa_dir, "weights",
	                                  "%s_%s-fov_wgt" % ( subj_id, id )
	                                )

	mvpa[ "wgt_ecc" ] = os.path.join( mvpa_dir, "weights",
	                                  "%s_%s-ecc_wgt" % ( subj_id, id )
	                                )

	paths[ "mvpa" ] = mvpa

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

	paths = _get_ana_paths( conf, paths )

	paths = _get_roi_paths( conf, paths )

	paths = _get_loc_paths( conf, paths )

	paths = _get_task_paths( conf, paths )

	paths = _get_mvpa_paths( conf, paths )

	return paths


def get_group_paths( conf ):
	"""Get the path structure for the group analysis"""

	exp_id = conf[ "exp" ][ "id" ]

	paths = {}

	paths[ "study" ] = _get_study_paths()

	paths[ "base_dir" ] = os.path.join( paths[ "study" ][ "base_dir" ],
	                                    "group_data"
	                                  )

	paths[ "uni" ] = os.path.join( paths[ "base_dir" ],
	                               "%s_group-uni.txt" % exp_id
	                             )

	return paths
