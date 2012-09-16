"""
Set of routines to run MVPA on the natural scenes
apertures fMRI experiment.
"""

from __future__ import division

import numpy as np


def task_xtr( paths, conf ):
	"""Extracts the details on task performance"""

	task = []

	# different items are different 'hand's
	keys = [ { "r" : 0,
	           "g" : 1,
	           "y" : 2,
	           "b" : 3
	         },
	         { "r" : 3,
	           "g" : 2,
	           "y" : 1,
	           "b" : 0
	         }
	       ]

	# order was ordered for subject s1000, and the order was not stored in the
	# task file
	s1000_hand_order = [ 1, 0, 0, 1, 1, 0, 0, 1, 1, 0 ]

	for i_run in xrange( conf[ "subj" ][ "n_exp_runs" ] ):

		task_file = "%s%d.npy" % ( paths[ "log" ][ "task_base" ], i_run + 1 )

		run_tasks = np.load( task_file )

		seq_file = "%s%d.npy" % ( paths[ "log" ][ "seq_base" ], i_run + 1 )

		run_seq = np.load( seq_file )

		# potential times that a task could have been cued
		evt_times = run_seq[ :, 0 ] + conf[ "task" ][ "evt_on_s" ]

		for curr_task in run_tasks:

			if conf[ "subj" ][ "subj_id" ] == "s1000":
				( resp, time_s ) = curr_task
				hand_flag = s1000_hand_order[ i_run ]
			else:
				( resp, time_s, hand_flag ) = curr_task

			# s1032 had it backwards
			if conf[ "subj" ][ "subj_id" ] == "s1032":
				hand_flag = np.logical_not( hand_flag )

			# find the last event where the response time is greater than the
			# presentation onset time
			i_evt = np.where( time_s > evt_times )[ 0 ][ -1 ]

			( cond, img_id_L, img_id_R ) = run_seq[ i_evt, [ 2, 5, 6 ] ]

			task.append( [ i_run,
			               cond,
			               keys[ hand_flag ][ resp ],
			               hand_flag,
			               img_id_L,
			               img_id_R
			             ]
			           )

	np.save( paths[ "task" ][ "resp" ], task )


def task_analysis( paths, conf ):
	"""Analyses the distribution of task responses"""

	resp = np.load( paths[ "task" ][ "resp" ] )

	hist = []

	for i_cond in xrange( 2 ):

		i_cond_resp = ( resp[ :, 1 ] == i_cond )

		# pull out the responses for this condition
		cond_resp = resp[ i_cond_resp, 2 ]

		hist.append( [ sum( cond_resp == button )
		               for button in [ 0, 1, 2, 3 ]
		             ]
		           )

	# button x cond
	hist = np.array( hist ).T

	# sum over buttons
	cond_sum = np.sum( hist, axis = 0 ).astype( "float" )

	# for some reason, /= does not obey the float division rules here, so do it
	# explicitly
	hist = hist / cond_sum

	np.save( paths[ "task" ][ "hist" ], hist )
