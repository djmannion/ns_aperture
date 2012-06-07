"""Natural scenes through apertures fMRI experiment.
"""

from __future__ import division

import os

import numpy as np

import psychopy.visual, psychopy.filters, psychopy.misc, psychopy.event
import psychopy.core, psychopy.log

import psychopy.logging
psychopy.logging.console.setLevel( psychopy.logging.CRITICAL )


import ns_aperture.config
import stimuli.psychopy_ext, stimuli.utils

__author__ = "Damien Mannion"
__license__ = "GPL"
__version__ = "1"
__maintainer__ = "Damien Mannion"
__email__ = "dmannion@umn.edu"
__status__ = "Collection"


def run( subj_id, run_num, order ):
	"""Execute a run of the natural scenes aperture fMRI experiment.

	Parameters
	----------
	subj_id : string
		Subject ID, in the Olman lab system.
	run_num : integer, > 0
		Run number.
	order : string, { "AB", "BA" }
		Block ordering.

	Returns
	-------
	status : integer, { 0, 1 }
		Returns 0 if run completed, 1 if aborted.

	"""

	# get the paths that will hold the run sequence and task information,
	# checking that they don't exist already
	( seq_path, task_path ) = get_log_paths( subj_id, run_num )

	# load the experiment configuration
	conf = ns_aperture.config.get_conf()

	# get the event sequence info
	seq = get_seq( conf, order )

	# get the dictionary that tells us what the columns in the sequence mean
	seq_ind = get_seq_ind()

	# initialise the behavioural task details
	task = init_task( conf )
	task_resp = []

	# initialise the display window
	win = psychopy.visual.Window( ( 1024, 768 ),
	                              monitor = conf[ "acq" ][ "monitor_name" ],
	                              fullscr = True,
	                              allowGUI = False
	                            )

	stat_txt = psychopy.visual.TextStim( win = win,
	                                     text = "Loading stimuli...",
	                                     height = 16,
	                                     units = "pix"
	                                   )

	# using random_integers because it is [ 0, 1 ] rather than [ 0, 1 )
	task_hand = np.random.random_integers( 0, 1, 1 )[ 0 ]

	if task_hand == 0:
		hand_instruct = "Same scene Left, different scene R"
	else:
		hand_instruct = "Same scene Right, different scene L"

	hand_txt = psychopy.visual.TextStim( win = win,
	                                     text = hand_instruct,
	                                     height = 16,
	                                     units = "pix",
	                                     pos = ( 0, -40 )
	                                   )

	try:
		fixation = init_fixation( conf, win )
		stim = init_stim( conf, win, stat_txt, hand_txt )
	except:
		win.close()
		raise

	# initialise the clock that keeps track of how long the run has been going
	run_clock = psychopy.core.Clock()

	# set the keys
	quit_key = 'q'
	trigger_key = 't'

	hand_txt.draw()

	# wait for the trigger
	fixation.draw()

	win.flip()

	k = psychopy.event.waitKeys( keyList = [ quit_key, trigger_key ] )

	# start the clock
	run_clock.reset()

	# quit, if the user wanted to
	if quit_key in k:
		print "User aborted"
		win.close()
		return 1

	run_time = run_clock.getTime()

	trig_count = 1

	# keep looping until the time has elapsed
	while run_time < conf[ "exp" ][ "run_len_s" ]:

		i_evt = np.where( run_time > seq[ :, seq_ind[ "time_s" ] ] )[ 0 ][ -1 ]

		evt_time = run_time - seq[ i_evt, seq_ind[ "time_s" ] ]

		evt_stims = [ stim[ int( seq[ i_evt, seq_ind[ "img_i_L" ] ] ) ][ 0 ],
		              stim[ int( seq[ i_evt, seq_ind[ "img_i_R" ] ] ) ][ 1 ]
		            ]

		# if the event time is less than how long we want the stimulus up for, draw
		# the stimuli
		if evt_time < conf[ "exp" ][ "evt_stim_s" ]:
			_ = [ evt_stim.draw() for evt_stim in evt_stims ]

		# set the fixation colour based on the event time, gated by whether it is a
		# task event or not
		if ( ( conf[ "task" ][ "evt_on_s" ] <
		       evt_time <
		       conf[ "task" ][ "evt_off_s" ]
		     ) and
		     ( i_evt in task )
		   ):
			fixation.setFillColor( conf[ "stim" ][ "fix_col_act" ] )
		else:
			fixation.setFillColor( conf[ "stim" ][ "fix_col_inact" ] )

		# draw the fixation
		fixation.draw()

		# draw to the screen
		win.flip()

		# get any responses
		keys = psychopy.event.getKeys( timeStamped = run_clock )

		for ( key, timestamp ) in keys:

			if key == quit_key:
				print "User abort"
				win.close()
				return 1

			elif key == trigger_key:
				trig_count += 1

			else:
				task_resp.extend( [ ( key, timestamp, task_hand ) ] )

		run_time = run_clock.getTime()

	# all done, time for cleanup
	# first, close the window
	win.close()

	# now save the sequence
	np.save( seq_path, seq )

	# convert the response list into an np array for saving
	task_resp_np = np.array( task_resp,
	                         dtype = [ ( "key", "S10" ),
	                                   ( "time", float )
	                                   ( "hand", int )
	                                 ]
	                       )
	# and save
	np.save( task_path, task_resp_np )

	return 0


def init_fixation( conf, win ):
	"""Initialises the fixation display.

	Parameters
	----------
	conf : dict
		Experiment configuration, as output from ``ns_apertures.config.get_conf()``
	win : Window object
		PsychoPy window.

	Returns
	-------
	fixation : Circle object
		PsychoPy stimulus object for fixation.

	"""

	fix_rad_pix = psychopy.misc.deg2pix( conf[ "stim" ][ "fix_rad_deg" ],
	                                     win.monitor
	                                   )

	fixation = psychopy.visual.Circle( win = win,
	                                   radius = fix_rad_pix,
	                                   units = "pix",
	                                   fillColor = conf[ "stim" ][ "fix_col_inact" ],
	                                   lineWidth = 1
	                                 )

	return fixation


def init_stim( conf, win, stat_txt, hand_txt ):
	"""Initialises the stimuli for the experiment.

	Parameters
	----------
	conf : dict
		Experiment configuration, as output from ``ns_apertures.config.get_conf()``
	win : Window object
		PsychoPy window.
	stat_txt : TextStim object
		PsychoPy text object

	Returns
	-------
	stim : list of PatchStim-s
		Outer list is for a given image, inner list is for the image regions.

	"""

	stim_conf = conf[ "stim" ]

	im_patch_d = np.round( stim_conf[ "patch_diam_deg" ] *
	                       ( 1.0 / stim_conf[ "im_deg_pp" ] )
	                     )

	mask_diam = stimuli.utils.nearest_power_of_two( im_patch_d )

	# make the aperture
	mask = psychopy.filters.makeMask( mask_diam,
	                                  radius = float( im_patch_d / mask_diam ),
	                                  shape = "circle",
	                                  fringeWidth = 0.0
	                                )

	# stimuli in the configuration are in the space of the images
	# now, we need to convert them to the space of our display

	# degrees per pixel for this monitor
	m_dpp = psychopy.misc.pix2deg( 1, win.monitor )

	# scale factor relative to the dpp of the image space
	scale_fac = stim_conf[ "im_deg_pp" ] / m_dpp

	patch_size = mask.shape[ 0 ] * scale_fac

	patch_pos = [ psychopy.misc.deg2pix( stim_conf[ "patch_ecc_deg" ], win.monitor ) *
	              offset
	              for offset in ( -1, +1 )
	            ]

	n = stim_conf[ "img_ids" ].shape[ 0 ]

	stim = []

	for ( i_im, im_id ) in enumerate( stim_conf[ "img_ids" ] ):

		stat_txt.setText( "Loading stimuli (%d/%d)" % ( i_im + 1, n ) )
		stat_txt.draw()

		hand_txt.draw()

		win.flip()

		im_stim = []

		for ( p_rect, p_pos ) in zip( stim_conf[ "patch_rect" ], patch_pos ):

			img = stimuli.utils.read_van_hateren( im_id,
			                                      path = stim_conf[ "db_path" ],
			                                      pad = True,
			                                      flip = True,
			                                      rect = p_rect,
			                                      scale = stim_conf[ "scale_mode" ]
			                                    )

			tex = psychopy.visual.PatchStim( win = win,
			                                 tex = img,
			                                 size = patch_size,
			                                 units = "pix",
			                                 mask = mask,
			                                 pos = ( p_pos, 0 )
			                               )

			im_stim.extend( [ tex ] )

		stim.append( im_stim )

	return stim


def get_seq_ind():
	"""Defines a lookup table for the columns in the experiment sequence array.

	Returns
	-------
	seq_ind : dict, containing items (all ints):
		time_s : event onset time, in seconds
		block_num : run block number
		block_type : whether an 'A' or 'B' block
		img_i_L : subset image index, for L stim
		img_i_R : subset image index, for R stim
		img_id_L : van Hateren image index, for L stim
		img_id_R : van Hateren image index, for R stim

	"""

	seq_ind = { "time_s" : 0,
	            "block_num" : 1,
	            "block_type" : 2,
	            "img_i_L" : 3,
	            "img_i_R" : 4,
	            "img_id_L" : 5,
	            "img_id_R" : 6
	          }

	return seq_ind


def get_seq( conf, order ):
	"""Get a sequence of events that consitute a run.

	Returns
	-------
	seq : ( evt, parameter ) numpy array
		Details for each event. The columns are given by ``get_seq_ind``

	"""

	seq_ind = get_seq_ind()

	# init the empty sequence
	seq = np.empty( ( conf[ "exp" ][ "n_evt_per_run" ],
	                  len( seq_ind )
	                )
	              )
	seq.fill( np.NAN )

	# repeat the image indices three times; first for the A blocks, second and
	# third for the B blocks
	img_i_set = np.tile( np.arange( conf[ "stim" ][ "img_ids" ].shape[ 0 ] ),
	                     ( 3, 1 )
	                   )

	# shuffle the indices
	map( np.random.shuffle, img_i_set )

	# loop through each event
	for i_evt in xrange( conf[ "exp" ][ "n_evt_per_run" ] ):

		# onset time
		time_s = i_evt * conf[ "exp" ][ "evt_len_s" ]

		# block number (one-indexed)
		block_num = np.floor( i_evt / conf[ "exp" ][ "n_evt_per_block" ] ) + 1

		# if the block is odd
		if np.mod( block_num, 2 ) == 1:

			if order == "AB":
				block_type = 0
			else:
				block_type = 1

		# if the block is even
		else:

			if order == "AB":
				block_type = 1
			else:
				block_type = 0

		seq[ i_evt, seq_ind[ "time_s" ] ] = time_s
		seq[ i_evt, seq_ind[ "block_num" ] ] = block_num
		seq[ i_evt, seq_ind[ "block_type" ] ] = block_type

	# find all the A events, and set them to the first image set for both the
	# left and the right patches
	a_evts = seq[ :, seq_ind[ "block_type" ] ] == 0
	seq[ a_evts, seq_ind[ "img_i_L" ] ] = img_i_set[ 0, : ].copy()
	seq[ a_evts, seq_ind[ "img_i_R" ] ] = img_i_set[ 0, : ].copy()

	# find all the B events, and set the L and R image sets to be different from
	# A and different from one another
	b_evts = seq[ :, seq_ind[ "block_type" ] ] == 1
	seq[ b_evts, seq_ind[ "img_i_L" ] ] = img_i_set[ 1, : ].copy()
	seq[ b_evts, seq_ind[ "img_i_R" ] ] = img_i_set[ 2, : ].copy()

	# get the list of image identifiers
	id_db = conf[ "stim" ][ "img_ids" ]

	# store the image identifier corresponding to each image index
	seq[ :, seq_ind[ "img_id_L" ] ] = id_db[ seq[ :, seq_ind[ "img_i_L" ] ].astype( "int" ) ]
	seq[ :, seq_ind[ "img_id_R" ] ] = id_db[ seq[ :, seq_ind[ "img_i_R" ] ].astype( "int" ) ]


	# some sanity checks
	# * filled up the array
	assert( np.sum( np.isnan( seq ) ) == 0 )

	# * same number of A and B events
	assert( ( np.sum( seq[ :, seq_ind[ "block_type" ] ] == 0 ) ==
	          np.sum( seq[ :, seq_ind[ "block_type" ] ] == 1 )
	        )
	      )

	# * left image set contains all images
	assert( ( np.unique( seq[ :, seq_ind[ "img_i_L" ] ] ).shape[ 0 ] ==
	          conf[ "stim" ][ "img_ids" ].shape[ 0 ]
	        )
	      )

	# * right image set contains all images
	assert( ( np.unique( seq[ :, seq_ind[ "img_i_R" ] ] ).shape[ 0 ] ==
	          conf[ "stim" ][ "img_ids" ].shape[ 0 ]
	        )
	      )

	return seq


def get_log_paths( subj_id, run_num ):
	"""Return the paths to run log files.

	Parameters
	----------
	subj_id : string
		Subject ID
	run_num : int (>0)
		Run number

	Returns
	-------
	seq_path : string
		Path to save the run image sequence.
	task_path : string
		Path to save the task and response sequence.

	Notes
	-----
	* If either of the files exist already, then an exception is thrown

	"""

	seq_path = os.path.join( "../logs",
	                         ( "%s_ns_aperture_fmri_seq_%d.npy" %
	                           ( subj_id, run_num )
	                         )
	                       )

	if os.path.exists( seq_path ):
		raise IOError( "".join( [ "Path ", seq_path, " already exists" ] ) )

	task_path = os.path.join( "../logs",
	                          ( "%s_ns_aperture_fmri_task_%d.npy" %
	                            ( subj_id, run_num )
	                          )
	                        )

	if os.path.exists( task_path ):
		raise IOError( "".join( [ "Path ", task_path, " already exists" ] ) )

	return ( seq_path, task_path )


def init_task( conf ):
	"""Initialises the task timing.

	Returns
	-------
	task_evts : numpy array of integers
		Indices for events where subjects should be asked to perform the task.

	"""

	# generate a set of random event deltas (time between events, in units of
	# events) from a geometric distribution
	# the size argument is arbitrary, just needs to be enough to have plenty of
	# headroom
	evt_delta = np.random.geometric( p = conf[ "task" ][ "p" ],
	                                 size = 500
	                               )

	# first task event occurs at the first delta
	task_evts = [ evt_delta[ 0 ] ]

	# recursively define future events as additional to the time of the previous
	# event
	task_evts.extend( ( task_evts[ -1 ] + delta )
	                  for delta in evt_delta
	                )

	# easier to handle if it is an array
	task_evts = np.array( task_evts )

	# restrict the events to those falling within the duration of a run
	task_evts = task_evts[ task_evts <= conf[ "exp" ][ "n_evt_per_run" ] ]

	return task_evts
