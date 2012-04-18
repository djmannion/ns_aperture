
import numpy as np

def run( subj_id, run_no, order ):
	"""
	"""

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
	task_evts.extend( ( task_events[ -1 ] + delta )
	                  for delta in evt_delta
	                )

	# easier to handle if it is an array
	task_evts = np.array( task_evts )

	# restrict the events to those falling within the duration of a run
	task_evts = task_evts[ task_evts <= conf[ "exp" ][ "n_events_per_run" ] ]

	return task_evts
