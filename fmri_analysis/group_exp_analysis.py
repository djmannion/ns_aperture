"""
Set of routines to analyse single-subject fMRI data for the natural scenes
apertures fMRI experiment.
"""

from __future__ import division

import os, os.path

import numpy as np

import fmri_tools.utils

import ns_aperture.fmri.exp


def agg( grp_paths, conf ):
	"""Aggregate the univariate analyses across subjects"""

	n_subj = len( conf[ "all_subj" ] )

	n_rois = len( conf[ "ana" ][ "rois" ] ) * len( conf[ "ana" ][ "i_parc" ] )

	uni_data = np.empty( ( n_subj, n_rois ) )
	uni_data.fill( np.NAN )

	for ( i_subj, subj_id ) in enumerate( conf[ "all_subj" ] ):

		subj_conf = ns_aperture.config.get_conf( subj_id )
		subj_paths = ns_aperture.fmri_analysis.paths.get_subj_paths( subj_conf )

		i_roi = 0

		for ( roi_name, _ ) in conf[ "ana" ][ "rois" ]:

			roi_psc = np.loadtxt( "%s_%s.txt" % ( subj_paths[ "rois" ][ "parc_psc" ],
			                                      roi_name
			                                    )
			                    )

			uni_data[ i_subj, i_roi:( i_roi + 2 ) ] = roi_psc

			i_roi += 2

	assert( np.sum( np.isnan( uni_data ) ) == 0 )

	np.savetxt( grp_paths[ "uni" ], uni_data )
