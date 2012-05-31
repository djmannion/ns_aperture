"""Natural scenes through apertures fMRI experiment localiser.
"""

from __future__ import division

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


def run( order ):
	"""Execute a localiser run for the natural scenes aperture fMRI experiment.

	Parameters
	----------
	order : string, { "AB", "BA" }
		Block ordering.

	Returns
	-------
	status : integer, { 0, 1 }
		Returns 0 if run completed, 1 if aborted.

	"""

	# load the experiment configuration
	conf = ns_aperture.config.get_conf()

	# get the event sequence info
	seq = get_seq( conf, order )

	# get the dictionary that tells us what the columns in the sequence mean
	seq_ind = get_seq_ind()

	# initialise the display window
	win = psychopy.visual.Window( ( 1024, 768 ),
	                              monitor = conf[ "acq" ][ "monitor_name" ],
	                              fullscr = True,
	                              allowGUI = False
	                            )

	try:
		fixation = init_fixation( conf, win )
		stim = init_stim( conf, win )
	except:
		win.close()
		raise

	# initialise the clock that keeps track of how long the run has been going
	run_clock = psychopy.core.Clock()

	# set the keys
	quit_key = 'q'
	trigger_key = 't'

	fix_text = psychopy.visual.TextStim( win = win,
	                                     text = "Awaiting trigger ...",
	                                     height = 16,
	                                     units = "pix",
	                                     bold = False,
	                                     pos = ( 0, -40 )
	                                   )

	fix_text.draw()

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

	# keep looping until the time has elapsed
	while run_time < conf[ "exp" ][ "run_len_s" ]:

		i_evt = np.where( run_time > seq[ :, seq_ind[ "time_s" ] ] )[ 0 ]

		if len( i_evt ) > 0:

			i_evt = i_evt[ -1 ]

			i_off = int( seq[ i_evt, seq_ind[ "block_type" ] ] ) - 1
			i_phase = int( seq[ i_evt, seq_ind[ "i_phase" ] ] )

			if i_off >= 0:

				evt_stim = stim[ i_off ][ i_phase ]

				evt_stim.draw()

		# draw the fixation
		fixation.draw()

		# draw to the screen
		win.flip()

		# get any responses
		keys = psychopy.event.getKeys( timeStamped = run_clock )

		for ( key, _ ) in keys:

			if key == quit_key:
				print "User abort"
				win.close()
				return 1

		run_time = run_clock.getTime()

	# all done, time for cleanup
	# first, close the window
	win.close()

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


def init_stim( conf, win ):
	"""Initialises the stimuli for the experiment.

	Parameters
	----------
	conf : dict
		Experiment configuration, as output from ``ns_apertures.config.get_conf()``
	win : Window object
		PsychoPy window.

	Returns
	-------
	stim : list of PatchStim-s
		Outer list is for a given image, inner list is for the image regions.

	"""

	stim_conf = conf[ "stim" ]

	patch_diam_pix = psychopy.misc.deg2pix( stim_conf[ "patch_diam_deg" ],
	                                        win.monitor
	                                      )

	tex_diam_pix = stimuli.utils.nearest_power_of_two( patch_diam_pix )

	tex_cycles = ( psychopy.misc.pix2deg( tex_diam_pix, win.monitor ) *
	               stim_conf[ "loc_sf_cpd" ]
	             )

	cheqs = [ psychopy.filters.makeGrating( tex_diam_pix,
	                                        ori = ori,
	                                        cycles = tex_cycles,
	                                        gratType = "sqr"
	                                      )
	         for ori in ( 0, 90 )
	        ]

	cheq_pA = cheqs[ 0 ] * cheqs[ 1 ]

	cheq_pB = ( cheqs[ 0 ] * -1 ) * cheqs[ 1 ]

	# make the aperture
	mask = psychopy.filters.makeMask( tex_diam_pix,
	                                  radius = float( patch_diam_pix / tex_diam_pix ),
	                                  shape = "circle",
	                                  fringeWidth = 0.0
	                                )

	patch_pos = [ psychopy.misc.deg2pix( stim_conf[ "patch_ecc_deg" ], win.monitor ) *
	              offset
	              for offset in ( -1, +1 )
	            ]


	stim = [ [ psychopy.visual.PatchStim( win = win,
	                                      tex = cheq_pA.copy(),
	                                      size = tex_diam_pix,
	                                      units = "pix",
	                                      mask = mask,
	                                      pos = ( p_pos, 0 )
	                                    ),
	           psychopy.visual.PatchStim( win = win,
	                                      tex = cheq_pB.copy(),
	                                      size = tex_diam_pix,
	                                      units = "pix",
	                                      mask = mask,
	                                      pos = ( p_pos, 0 )
	                                    )
	        ]
	        for p_pos in patch_pos
	      ]

	return stim


def get_seq_ind():
	"""Defines a lookup table for the columns in the experiment sequence array.

	Returns
	-------
	seq_ind : dict, containing items (all ints):
		time_s : event onset time, in seconds
		block_num : run block number
		block_type : whether an 'A' or 'B' block
		contrast_L : contrast of the left stimulus
		contrast_R : contrast of the right stimulus
		i_phase : index of the phase of the stimulus

	"""

	seq_ind = { "time_s" : 0,
	            "block_num" : 1,
	            "block_type" : 2,
	            "contrast_L" : 3,
	            "contrast_R" : 4,
	            "i_phase" : 5,
	          }

	return seq_ind

def get_block_seq( order, n_blocks ):
	"""
	"""

	if order == "AB":
		mini_block = np.array( [ 1, 2, 0 ] )
	elif order == "BA":
		mini_block = np.array( [ 2, 1, 0 ] )

	n_reps = int( ( n_blocks - 1 ) / len( mini_block ) )

	# hard-coded number of blocks
	blk_seq = np.concatenate( ( [ 0 ], np.tile( mini_block, n_reps ) ) )

	return blk_seq


def get_seq( conf, order ):
	"""Get a sequence of events that consitute a run.

	Returns
	-------
	seq : ( evt, parameter ) numpy array
		Details for each event. The columns are given by ``get_seq_ind``

	"""

	block_seq = get_block_seq( order, conf[ "exp" ][ "loc_n_blocks" ] )

	seq_ind = get_seq_ind()

	n_evt = conf[ "stim" ][ "loc_rev_rate_hz" ] * conf[ "exp" ][ "run_len_s" ]

	# init the empty sequence
	seq = np.empty( ( n_evt,
	                  len( seq_ind )
	                )
	              )
	seq.fill( np.NAN )

	i_phase = 0

	for i_evt in xrange( int( n_evt ) ):

		time_s = i_evt * ( 1.0 / conf[ "stim" ][ "loc_rev_rate_hz" ] )

		time_s += conf[ "exp" ][ "loc_pre_len_s" ]

		i_block = np.floor( i_evt *
		                    ( 1.0 / conf[ "stim" ][ "loc_rev_rate_hz" ] ) /
		                    conf[ "exp" ][ "block_len_s" ]
		                  )

		block_type = block_seq[ i_block ]

		if block_type == 0:
			contrast_L = 0.0
			contrast_R = 0.0

		elif block_type == 1:
			contrast_L = 1.0
			contrast_R = 0.0

		elif block_type == 2:
			contrast_L = 0.0
			contrast_R = 1.0

		i_phase = int( not( i_phase ) )

		seq[ i_evt, seq_ind[ "time_s" ] ] = time_s
		seq[ i_evt, seq_ind[ "block_num" ] ] = i_block + 1
		seq[ i_evt, seq_ind[ "block_type" ] ] = block_type
		seq[ i_evt, seq_ind[ "contrast_L" ] ] = contrast_L
		seq[ i_evt, seq_ind[ "contrast_R" ] ] = contrast_R
		seq[ i_evt, seq_ind[ "i_phase" ] ] = i_phase


	# some sanity checks
	# * filled up the array
	assert( np.sum( np.isnan( seq ) ) == 0 )

	return seq
