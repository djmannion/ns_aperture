

import numpy as np


class ConfigContainer( object ):

	pass


def get_conf( subj_id = None ):
	"""Overall experiment configuration"""

	conf = ConfigContainer()

	conf.stim = _get_stim_conf()
	conf.acq = _get_acq_conf()
	conf.exp = _get_exp_conf()
	conf.ana = _get_ana_conf( conf )
	conf.task = _get_task_conf()
	conf.all_subj = _get_subj_conf( None )

	if subj_id is not None:
		conf = _apply_subj_specific_fixes( conf, subj_id )
		conf.subj = _get_subj_conf( subj_id )

	return conf


def _get_stim_conf():
	"""Get the stimulus configuration"""

	stim_conf = ConfigContainer()

	# path to the image database (the image files)
	stim_conf.db_path = "../im_db"

	# ( n_rows, n_cols ) dimensions of each (full) image
	stim_conf.im_dim = ( 1024, 1536 )

	# degrees per pixel in the images as acquired (ie native)
	stim_conf.im_deg_pp = 1.0 / 60.0

	# diameter, in degrees visual angle, of each extracted image patch
	stim_conf.patch_diam_deg = 4.0

	# eccenticity, in degrees visual angle, where each extracted patch will be
	# shown
	stim_conf.patch_ecc_deg = 3.0

	patch_diam_pix = np.round( stim_conf.patch_diam_deg *
	                           ( 1.0 / stim_conf.im_deg_pp )
	                         )

	patch_r_pix = patch_diam_pix / 2

	patch_ecc_pix = np.round( stim_conf.patch_ecc_deg *
	                          ( 1.0 / stim_conf.im_deg_pp )
	                        )

	patch_centre_pix = [ ( ( stim_conf.im_dim[ 0 ] / 2 ),
	                       ( stim_conf.im_dim[ 1 ] / 2 +
	                         offset * patch_ecc_pix
	                       )
	                     )
	                     for offset in ( -1, +1 )
	                   ]

	# coordinates of each patch, in the original image dimensions. Arranged as
	# ( patch 1, patch 2 ) -> ( top left, bottom right ) -> ( i_row, i_col )
	stim_conf.patch_rect = [ ( ( c_pix[ 0 ] - patch_r_pix, # top
	                             c_pix[ 1 ] - patch_r_pix  # left
	                           ),
	                           ( c_pix[ 0 ] + patch_r_pix, # bottom
	                             c_pix[ 1 ] + patch_r_pix  # right
	                           )
	                         )
	                         for c_pix in patch_centre_pix
	                       ]

	# how to scale the images.
	stim_conf.scale_mode = "mean"

	# hardcode the desired image ids
	stim_conf.img_ids = np.array( (    6,   12,  504,   58,   72,   90,
	                                  98,  110,  181,  230,  285,  342,
	                                 400,  460,  474,  482,  485,  561,
	                                 567,  639,  664,  727,  730,  737,
	                                 752,  756,  821, 1068, 1083, 1087,
	                                1117, 1166, 1234, 1270, 1398, 1402,
	                                1435, 1448, 1466, 1467, 1471, 1472,
	                                1489, 1501, 1520, 1523, 1526, 1531,
	                                1533, 1537, 1543, 1544, 1559, 1565,
	                                1566, 1568, 1570, 1574, 1575, 1578,
	                                1579, 1580, 1583, 1584, 1601, 1611,
	                                1614, 1618, 1626, 1629, 1632, 1639,
	                                1647, 1654, 1656, 1657, 1663, 1676,
	                                1679, 1933, 1950, 1957, 2066, 2260,
	                                2329, 2339, 2360, 2429, 2745, 2889,
	                                2994, 2997, 3024, 3035, 3057, 3061,
	                                3065, 3086, 3099, 3194, 3261, 3266,
	                                3276, 3287, 3297, 3302, 3325, 3336
	                             )
	                           )


	# radius of the fixation circle, in degrees visual angle.
	stim_conf.fix_rad_deg = 0.075

	# colour of inner fixation circle when task is inactive.
	stim_conf.fix_col_inact = ( -1, -1, -1 )

	# colour of inner fixation circle when task is active.
	stim_conf.fix_col_act = ( -1, 1, -1 )

	# spatial frequency of the localiser stimulus, in cycles per degree
	stim_conf.loc_sf_cpd = 2.0

	# reversal rate of the localiser stimulus, in hertz.
	stim_conf.loc_rev_rate_hz = 2

	return stim_conf


def _get_exp_conf():
	"""Gets the experiment configuration"""

	exp_conf = ConfigContainer()

	exp_conf.id = "ns_aperture"
	exp_conf.id_loc = "ns_aperture_loc"

	# main experiment info
	exp_conf.n_runs = 10
	exp_conf.n_blocks = 18
	exp_conf.block_len_s = 16.0
	exp_conf.run_len_s = exp_conf.n_blocks * exp_conf.block_len_s


	exp_conf.n_evt_per_block = 12
	exp_conf.n_evt_per_run = exp_conf.n_evt_per_block * exp_conf.n_blocks
	exp_conf.evt_len_s = exp_conf.block_len_s / exp_conf.n_evt_per_block

	exp_conf.evt_stim_s = 1.0

	# localiser info
	exp_conf.loc_n_runs = 2
	exp_conf.loc_n_blocks = 19
	exp_conf.loc_run_len_s = exp_conf.loc_n_blocks * exp_conf.block_len_s
	exp_conf.loc_pre_len_s = 6

	exp_conf.loc_run_full_len_s = exp_conf.loc_run_len_s + exp_conf.loc_pre_len_s

	return exp_conf


def _get_task_conf():
	"""Gets the task configuration"""

	task_conf = ConfigContainer()

	task_conf.p = 0.35

	task_conf.evt_on_s = 0.8
	task_conf.evt_off_s = 1.2

	return task_conf


def _get_acq_conf():
	"""Get the acquisition configuration"""

	acq_conf = ConfigContainer()

	acq_conf.monitor_name = "UMN_7T"

	acq_conf.tr_s = 2.0

	# how to reshape the data to be in +RAS convention
	# subsequent commands are relative to the data AFTER this operation
	# see docs for how to determine these
	acq_conf.reshape_to_RAS = ( "-x", "-z", "-y" )

	# phase encode direction, according to the data's internal axes
	acq_conf.ph_encode_dir = "z"

	acq_conf.ph_enc_afni = "SI"

	acq_conf.fugue_params = [ "--median" ]

	acq_conf.vox_size = "1.5"

	acq_conf.vol_base = 72

	# TE difference in the fieldmaps
	acq_conf.delta_te_ms = 1.02

	# corresponds to the echo spacing
	acq_conf.dwell_ms = 0.65 / 2.0

	return acq_conf


def _get_ana_conf( conf ):
	"""Gets the parameters for the fMRI analysis"""

	ana_conf = ConfigContainer()

	# cull the first and last block
	ana_conf.exp_pre_cull_s = conf.exp.block_len_s
	ana_conf.exp_post_cull_s = conf.exp.block_len_s

	# cull the first few volumes
	ana_conf.loc_run_start_s = conf.exp.loc_pre_len_s
	ana_conf.loc_run_dur_s = conf.exp.loc_run_len_s

	# p threshold for the localiser
	ana_conf.loc_p = 0.01

	ana_conf.poly_ord = "a"

	ana_conf.hrf_model = "SPMG1({t:0f})".format( t = conf.exp.block_len_s )

	ana_conf.hrf_corr_vol = 3

	ana_conf.smooth_fwhm = "2.5"
	ana_conf.smooth_sigma = "1.0"

	ana_conf.rep_subj_id = "s1021"

	ana_conf.p_height_thr = 0.01

	ana_conf.mvpa_filt_ord = "2"
	ana_conf.mvpa_blocks = conf.exp.n_blocks - 2
	ana_conf.block_len_vol = int( conf.exp.block_len_s / conf.acq.tr_s )
	ana_conf.i_start_hrf_corr = ana_conf.block_len_vol + ana_conf.hrf_corr_vol

	ana_conf.sl_r = 5

	ana_conf.n_common_nodes = { "lh" : 44697, "rh" : 39020 }

	return ana_conf



def _get_subj_conf( subj_id = None ):
	"""Gets the configuration info for subjects"""

	s1000 = ConfigContainer()

	s1000.subj_id = "s1000"
	s1000.acq_date = "20120601"
	s1000.n_runs = 12
	s1000.n_exp_runs = 10
	s1000.n_loc_runs = 2
	s1000.mot_base = 7
	s1000.exp_runs =  [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
	s1000.loc_runs = [ 11, 12 ]
	s1000.comments = ""
	s1000.mask_SI = 90
	s1000.vis_loc_date = "20111212"

	s1000.extra_al_params = [ "-parang", "1", "-6", "4",
	                          "-parang", "2", "11", "21",
	                          "-parang", "3", "23", "33",
	                          "-maxrot", "10",
	                          "-source_automask+2",
	                          "-nocmass"
	                        ]

	s1000.node_k = { "lh" : 130318,
	                 "rh" : 131151
	               }

	s1021 = ConfigContainer()

	s1021.subj_id = "s1021"
	s1021.acq_date = "20120608"
	s1021.n_runs = 12
	s1021.n_exp_runs = 10
	s1021.n_loc_runs = 2
	s1021.mot_base = 7
	s1021.exp_runs =  [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
	s1021.loc_runs = [ 11, 12 ]
	s1021.comments = ""
	s1021.mask_SI = 90
	s1021.vis_loc_date = "20111212"

	s1021.extra_al_params = [ "-parang", "1", "-9", "1",
	                          "-parang", "2", "9", "19",
	                          "-parang", "3", "40", "54",
	                          "-maxrot", "10",
	                          "-source_automask+2",
	                          "-nocmass"
	                        ]

	s1021.node_k = { "lh" : 140847,
	                 "rh" : 141381
	               }

	s1011 = ConfigContainer()

	s1011.subj_id = "s1011"
	s1011.acq_date = "20120608"
	s1011.n_runs = 12
	s1011.n_exp_runs = 10
	s1011.n_loc_runs = 2
	s1011.mot_base = 7
	s1011.exp_runs =  [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
	s1011.loc_runs = [ 11, 12 ]
	s1011.comments = ""
	s1011.mask_SI = 90
	s1011.vis_loc_date = "20121219"

	s1011.extra_al_params = [ "-parang", "1", "-4", "6",
	                          "-parang", "2", "5", "15",
	                          "-parang", "3", "21", "31",
	                          "-maxrot", "10",
	                          "-source_automask+2",
	                          "-nocmass"
	                        ]

	s1011.node_k = { "lh" : 128434,
	                 "rh" : 128461
	               }

	s1032 = ConfigContainer()

	s1032.subj_id = "s1032"
	s1032.acq_date = "20120608"
	s1032.n_runs = 12
	s1032.n_exp_runs = 10
	s1032.n_loc_runs = 2
	s1032.mot_base = 7
	s1032.exp_runs =  [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
	s1032.loc_runs = [ 11, 12 ]
	s1032.comments = "Visible motion in last run (loc)"
	s1032.mask_SI = 90
	s1032.vis_loc_date = "20120203"

	s1032.extra_al_params = [ "-parini", "1", "13", #"18",
	                          "-parini", "2", "29", #"34",
	                          "-parini", "3", "22",# "27",
	                          "-maxrot", "10",
	                          "-source_automask+2",
	                          "-nocmass"
	                        ]

	s1032.node_k = { "lh" : 126078,
	                 "rh" : 126201
	               }

	s1008 = ConfigContainer()

	s1008.subj_id = "s1008"
	s1008.acq_date = "20120619"
	s1008.n_runs = 12
	s1008.n_exp_runs = 10
	s1008.n_loc_runs = 2
	s1008.mot_base = 7
	s1008.exp_runs =  [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
	s1008.loc_runs = [ 11, 12 ]
	s1008.comments = ""
	s1008.mask_SI = 90
	s1008.vis_loc_date = "20111214"

	s1008.extra_al_params = [ "-parang", "1", "0", "14",
	                          "-parang", "2", "11", "21",
	                          "-parang", "3", "0", "10",
	                          "-maxrot", "10",
	                          "-source_automask+2",
	                          "-nocmass"
	                        ]

	s1008.node_k = { "lh" : 140427,
	                 "rh" : 141898
	               }

	subj = ConfigContainer()

	subj.subj = { "s1000" : s1000,
	              "s1021" : s1021,
	              "s1011" : s1011,
	              "s1032" : s1032,
	              "s1008" : s1008
	            }

	if subj_id is None:
		return subj
	else:
		return subj.subj[ subj_id ]


def _apply_subj_specific_fixes( conf, subj_id ):
	"""Modify the experiment configuration for unexpected variations for
	particular subjects (hopefully these are rare"""


	if subj_id == "s1000":

		# for this subject (the first), the localiser run length didn't work
		# properly - it was ending at ``run_len_s`` rather than the
		# ``loc_run_full_len_s``, so we need to change the localiser duration
		# info.

		conf.exp.loc_n_blocks = 16
		conf.exp.loc_run_len_s = conf.exp.loc_n_blocks * conf.exp.block_len_s

		conf.exp.loc_run_full_len_s = conf.exp.loc_run_len_s + conf.exp.loc_pre_len_s

	return conf
