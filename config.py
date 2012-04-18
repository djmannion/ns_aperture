"""CONFIG
"""

from __future__ import division

import numpy as np

import fmri_tools.utils

def get_conf():
	"""
	"""

	conf = { "exp" : _get_exp_conf(),
	         "stim" : _get_stim_conf(),
	         "task" : _get_task_conf(),
	         "acq" : _get_acq_conf()
	       }

	return conf


def _get_stim_conf():
	"""Get the stimulus configuration.

	Specifies the stimulus configuration and parameters for the natural scenes
	aperture fMRI experiment.

	This shouldn't be called directly.

	Returns
	-------------
	db_path : string
		Path to the image database (the image files)
	db_sel_path : string
		Path to the numpy datafile specifying the image set
	im_dim : 2-item list
		( n_rows, n_cols ) dimensions of each (full) image
	im_deg_pp : scalar float
		Degrees per pixel in the images as acquired (ie native)
	patch_diam_deg : scalar float
		Diameter, in degrees visual angle, of each extracted image patch
	patch_ecc_deg : scalar float
		Eccenticity, in degrees visual angle, where each extracted patch will be
		shown
	patch_rect : 2-item list of 2-item list of 2-item list
		Coordinates of each patch, in the original image dimensions. Arranged as
		( patch 1, patch 2 ) -> ( top left, bottom right ) -> ( i_row, i_col )

	Notes
	-----
	* Returns are contained within a dictionary.

	"""

	stim_conf = {}

	stim_conf[ "db_path" ] = "../im_db"

	stim_conf[ "db_sel_path" ] = "../data/im_sel.npy"

	stim_conf[ "im_dim" ] = ( 1024, 1536 )
	stim_conf[ "im_deg_pp" ] = 1.0 / 60.0

	stim_conf[ "patch_diam_deg" ] = 4.0
	stim_conf[ "patch_ecc_deg" ] = 3.0

	patch_diam_pix = np.round( stim_conf[ "patch_diam_deg" ] *
	                           ( 1.0 / stim_conf[ "im_deg_pp" ] )
	                         )

	patch_r_pix = patch_diam_pix / 2

	patch_ecc_pix = np.round( stim_conf[ "patch_ecc_deg" ] *
	                          ( 1.0 / stim_conf[ "im_deg_pp" ] )
	                        )

	patch_centre_pix = [ ( ( stim_conf[ "im_dim" ][ 0 ] / 2 ),
	                       ( stim_conf[ "im_dim" ][ 1 ] / 2 +
	                         offset * patch_ecc_pix
	                       )
	                     )
	                     for offset in ( -1, +1 )
	                   ]

	stim_conf[ "patch_rect" ] = [ ( ( c_pix[ 0 ] - patch_r_pix, # top
	                                  c_pix[ 1 ] - patch_r_pix  # left
	                                ),
	                                ( c_pix[ 0 ] + patch_r_pix, # bottom
	                                  c_pix[ 1 ] + patch_r_pix  # right
	                                )
	                              )
	                              for c_pix in patch_centre_pix
	                            ]

	# hardcode the desired image ids
	stim_conf[ "img_ids" ] = np.array( (    6,   12,   47,   58,   72,   90,
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


	return stim_conf


def _get_exp_conf():
	"""Gets the experiment configuration.

	Specifies the experiment configuration and parameters for the natural scenes
	aperture fMRI experiment.

	This shouldn't be called directly.

	Returns
	-------
	n_runs : scalar integer
		Number of experiment runs per session.
	n_blocks : scalar integer
		Number of blocks per run.
	block_len_s : scalar float
		Length of each block, in seconds.
	run_dur_s : scalar float
		Length of each run, in seconds.
	n_evt_per_block : scalar integer
		Number of events per block.
	n_evt_per_run : scalar integer
		Number of events per run.
	evt_len_s : scalar float
		Length of each 'event' within a block, in seconds.

	Notes
	-----
	* Returns are contained within a dictionary.

	"""

	exp_conf = {}

	exp_conf[ "n_runs" ] = 10

	exp_conf[ "n_blocks" ] = 18

	exp_conf[ "block_len_s" ] = 16.0

	exp_conf[ "run_len_s" ] = exp_conf[ "n_blocks" ] * exp_conf[ "block_len_s" ]

	exp_conf[ "n_evt_per_block" ] = 12

	exp_conf[ "n_evt_per_run" ] = ( exp_conf[ "n_evt_per_block" ] *
	                                exp_conf[ "n_blocks" ]
	                              )

	exp_conf[ "evt_len_s" ] = exp_conf[ "block_len_s" ] / exp_conf[ "n_evt_per_block" ]

	return exp_conf


def _get_task_conf():
	"""Gets the task configuration.

	Specifies the configuration for the behavioural task in the natural scenes
	aperture fMRI experiment.

	This shouldn't be called directly.

	Returns
	-------
	p : scalar float
		Probability that a given event will be a task trial.

	Notes
	-----
	* Returns are contained within a dictionary.

	"""

	task_conf = {}

	task_conf[ "p" ] = 0.35

	return task_conf


def _get_acq_conf():
	"""Get the acquisition configuration.

	Specifies the acquisition configuration and parameters for the
	natural scenes aperture fMRI experiment.

	Returns
	-------
	monitor_name : string
		monitor configuration name.
	tr_s : float
		time-to-repetition (TR), in seconds.
	delta_te_ms : float
		echo time differences for the fieldmaps, in milliseconds.
	dwell_ms : float
		dwell time, in milliseconds.
	slice_order : array of int
		slice acqusition indices, where 0 is the first slice.

	Notes
	-----
	* Returns are contained within a dictionary.

	"""

	acq_conf = {}

	acq_conf[ "monitor_name" ] = "UMN_7T"

	acq_conf[ "tr_s" ] = 2.0

	acq_conf[ "delta_te_ms" ] = 1.02

	acq_conf[ "dwell_ms" ] = 0.325

	acq_conf[ "slice_order" ] = fmri_tools.utils.get_slice_order( 36 )

	return acq_conf


def get_subj_conf():
	"""Gets the configuration info for each subject.

	Returns
	-------
	subj_id : string
		Subject ID, in the Olman lab system.
	acq_date : string
		Acquisition date, in YYYYMMDD format.
	comments : string
		Any comments about the scanning session.

	Notes
	-----
	* Return values are within a dictionary, which is itself within a dictionary
	  indexed by subject ID.

	"""

#	s1000 = { "subj_id" : "s1000",
#	          "acq_date" : "20120106",
#	          "comments" : ""
#	        }

#	subj_conf = { "s1000" : s1000,
#	            }

	subj_conf = {}

	return subj_conf
