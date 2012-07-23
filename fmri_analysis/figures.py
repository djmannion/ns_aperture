"""
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.stats


def _set_defaults():
	"""Set some sane defaults for figures.
	"""

	params = { 'axes.labelsize': 8,
	           'font.family' : 'Arial',
	           'font.sans-serif' : 'Helvetica',
	           'text.fontsize': 12,
	           'legend.fontsize': 7,
	           'xtick.labelsize': 5,
	           'xtick.direction' : 'out',
	           'xtick.major.size' : 2,
	           'ytick.labelsize': 5,
	           'ytick.direction' : 'out',
	           'ytick.major.size' : 2
	         }
	
	plt.rcParams.update( params )

	plt.ioff()


def _cleanup_fig( ax ):
	"""Apply some standard commands to clean up the axes on figures.
	"""

	for loc, spine in ax.spines.iteritems():

		spine.set_linewidth( 0.5 )

		if loc in [ "left", "bottom" ]:
			spine.set_position( ( "outward", 5 ) )
		elif loc in [ "right", "top" ]:
			spine.set_color( "none" )
		else:
			raise ValueError( "Unknown spine location: %s" % loc )

	ax.xaxis.set_ticks_position( "bottom" )
	ax.yaxis.set_ticks_position( "left" )


def plot_vox_psc( paths, conf ):
	"""a"""

	_set_defaults()

	figs = [ [ "v1", "v2", "v3" ],
	         [ "hv4", "v3ab" ],
	         [ "lo1", "lo2", "loc" ],
	         [ "hmtp" ]
	       ]

	fig_lims = [ 2, 3, 2, 1.5 ]

	i_loc_F = 14

	crit_F = np.sqrt( scipy.stats.f.ppf( 0.999, 2, 292 ) )

	for ( i_fig, fig_rois ) in enumerate( figs ):

		fig = plt.figure()

		fig.set_size_inches( 8, 5, forward = False )

		gs = gridspec.GridSpec( 1, len( fig_rois ) )

		xlims = []
		ylims = []

		axes = []

		for ( i_subfig, fig_roi ) in enumerate( fig_rois ):

			ax = plt.subplot( gs[ i_subfig ] )

			loc_F = []
			exp_psc = []

			for hemi in [ "lh", "rh" ]:

				loc_file = "%s_%s_%s.txt" % ( paths[ "ana" ][ "loc_roi" ],
				                              fig_roi,
				                              hemi
				                            )

				loc_roi = np.loadtxt( loc_file )

				loc_F.append( loc_roi[ :, i_loc_F ].copy() )

				exp_file = "%s_%s_%s.txt" % ( paths[ "ana" ][ "exp_roi" ],
				                              fig_roi,
				                              hemi
				                            )

				exp_roi = np.loadtxt( exp_file )

				exp_psc.append( exp_roi.copy() )

			loc_F = np.concatenate( loc_F )
			exp_psc = np.concatenate( exp_psc )

			ax.scatter( np.sqrt( loc_F ),
			            exp_psc,
			            facecolor = "k",
			            edgecolor = "k",
			            alpha = 0.1
			          )

			ax.hold( True )

			ax.set_title( fig_roi.upper() )

			ylims.append( np.max( np.abs( ax.get_ylim() ) ) )

			axes.append( ax )

		ylim = np.max( ylims )

		for ax in axes:

			xlim = ax.get_xlim()

			ax.plot( xlim, [ 0, 0 ], "k--" )

			ax.plot( [ crit_F, crit_F ],
			         [ -ylim, ylim ],
			         "r--"
			       )

			ax.set_ylim( ( -ylim, ylim  ) )
			ax.set_xlim( xlim )

	plt.show()


def subj_cond_diff( paths, conf ):
	"""
	"""

	_set_defaults()

	fig = plt.figure()

	fig.set_size_inches( 4.8, 3, forward = True )

	ax = fig.gca()

	ax.hold( True )

	for ( i_roi, roi_name ) in enumerate( conf[ "ana" ][ "rois" ] ):

		diff = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "cond_diff" ],
		                                roi_name
		                              )
		              )

		ax.plot( ( i_roi, i_roi ), diff[ 1: ], "k" )

		ax.scatter( i_roi, diff[ 0 ],
		            facecolor = "k",
		            edgecolor = "w",
		            s = 35,
		            zorder = 100
		           )

	ax.plot( ax.get_xlim(), ( 0, 0 ), "k--" )

	_cleanup_fig( ax )

	ax.set_xlim( ( -0.5, len( conf[ "ana" ][ "rois" ] ) - 0.5 ) )

	ax.set_ylim( ( -0.5, 0.5 ) )
	ax.set_xticks( np.arange( len( conf[ "ana" ][ "rois" ] ) ) )
	ax.set_xticklabels( conf[ "ana" ][ "rois" ] )

	ax.set_ylabel( "Same - different scene (psc)" )
	ax.set_xlabel( "Visual area" )

	fig = plt.gcf()
	fig.set_size_inches( 4.8, 3 )

	fig.show()

	fig = plt.gcf()
	fig.set_size_inches( 4.8, 3 )

	plt.show()


def subj_resp_over_time( paths, conf ):
	"""
	"""

	_set_defaults()



	fig = plt.figure()

	fig.set_size_inches( 7, 5, forward = True )

	gs = gridspec.GridSpec( 2, 3 )

	for ( i_roi, roi_name ) in enumerate( conf[ "ana" ][ "rois" ] ):

		ax = plt.subplot( gs[ i_roi ] )

		ax.hold( True )

		blks = np.load( "%s-%s.npy" % ( paths[ "ana" ][ "block" ],
		                                roi_name
		                              )
		              )

		cond_a = blks[ blks[ :, 1 ] == 0, 0 ]
		cond_b = blks[ blks[ :, 1 ] == 1, 0 ]


		ax.plot( cond_a, "b" )
		ax.plot( cond_b, "r" )

		ax.plot( ax.get_xlim(), ( 0, 0 ), "k--", zorder = -1 )

		ax.text( 0.82,
		         0.9,
		         roi_name,
		         transform = ax.transAxes,
		         fontsize = 10 / 1.25
		       )

		_cleanup_fig( ax )

		ax.set_ylabel( "Response (psc)" )
		ax.set_xlabel( "Time (blocks)" )

	plt.subplots_adjust( left = 0.10,
	                     bottom = 0.10,
	                     right = 0.97,
	                     top = 0.90,
	                     wspace = 0.40,
	                     hspace = 0.34
	                   )

	plt.show()
