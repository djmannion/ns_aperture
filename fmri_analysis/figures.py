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

	left, width = 0.1, 0.8
	bottom, height = 0.1, 0.6
	bottom_h = left_h = left+height+0.02

	rect_scatter = [ left, bottom, width, height ]
	rect_hist = [ left, bottom_h, width, 0.23 ]

	i_loc_t = np.array( [ 2, 5 ] )

	crit_t = scipy.stats.t.ppf( 0.999, 292 )

	dot_cols = [ [ 0, 0, 1 ], [ 1, 0, 0 ] ]

	for ( i_roi, roi ) in enumerate( conf[ "ana" ][ "rois" ] ):

		( roi_name, _ ) = roi

		fig = plt.figure()
		fig.set_size_inches( 8, 5, forward = False )

		ax_scatter = plt.axes( rect_scatter)
		ax_hist = plt.axes( rect_hist )

		xlims = []
		ylims = []

		axes = []

		loc_psc = []
		exp_psc = []
		cols = []

		sig = []

		for hemi in [ "lh", "rh" ]:

			loc_stat_file = "%s_%s_%s.txt" % ( paths[ "ana" ][ "loc_roi_stat" ],
			                                   roi_name,
			                                   hemi
			                                 )

			loc_roi_stat = np.loadtxt( loc_stat_file )

			loc_psc_file = "%s_%s_%s.txt" % ( paths[ "ana" ][ "loc_roi_psc" ],
			                                  roi_name,
			                                  hemi
			                                )

			loc_roi_psc_all = np.loadtxt( loc_psc_file )

			i_psc = np.argmax( loc_roi_stat[ :, i_loc_t ], axis = 1 )
			t_max = np.max( loc_roi_stat[ :, i_loc_t ], axis = 1 )

			loc_roi_psc = np.empty( loc_roi_psc_all.shape[ 0 ] )

			for i_vox in xrange( loc_roi_psc_all.shape[ 0 ] ):

				loc_roi_psc[ i_vox ] = loc_roi_psc_all[ i_vox, i_psc[ i_vox ] ]

				if t_max[ i_vox ] >= crit_t:
					cols.append( dot_cols[ 0 ] )
					sig.append( 1 )
				else:
					cols.append( dot_cols[ 1 ] )
					sig.append( 0 )

			loc_psc.append( loc_roi_psc.copy() )

			exp_file = "%s_%s_%s.txt" % ( paths[ "ana" ][ "exp_roi_psc" ],
			                              roi_name,
			                              hemi
			                            )

			exp_roi = np.loadtxt( exp_file )

			exp_psc.append( exp_roi.copy() )

		loc_psc = np.concatenate( loc_psc )
		exp_psc = np.concatenate( exp_psc )

		sig = np.array( sig )

		ns_kde = scipy.stats.gaussian_kde( exp_psc[ sig == 0 ] )
		s_kde = scipy.stats.gaussian_kde( exp_psc[ sig == 1 ] )

		ax_scatter.scatter( exp_psc,
		                    loc_psc,
		                    facecolor = cols,
		                    edgecolor = cols,
		                    alpha = 0.1
		                  )

		xl = ax_scatter.get_xlim()

		xl_max = np.max( np.abs( xl ) )

		ax_scatter.set_xlim( ( -xl_max, xl_max ) )

		ax_scatter.set_xlabel( "Coherent scene beta (psc)" )
		ax_scatter.set_ylabel( "Localiser beta (psc)" )

#		ax_hist.hist( exp_psc[ sig == 0 ], orientation = "horizontal" )

		x = np.linspace( -xl_max, xl_max, 500 )

		ns_eval = ns_kde( x )
		s_eval = s_kde( x )

		ax_hist.plot( x, ns_eval, color=dot_cols[ 1 ] )
		ax_hist.hold( True )
		ax_hist.plot( x, s_eval, color=dot_cols[ 0 ] )

		ax_hist_ylim = ax_hist.get_ylim()

		ax_hist.plot( [ 0, 0 ], ax_hist_ylim, "k--" )

		ns_mean_psc = np.mean( exp_psc[ sig == 0 ] )
		s_mean_psc = np.mean( exp_psc[ sig == 1 ] )

		ax_hist.plot( [ ns_mean_psc, ns_mean_psc ],
		              ax_hist_ylim,
		              color = dot_cols[ 1 ],
		              linestyle = "--"
		            )

		ax_hist.plot( [ s_mean_psc, s_mean_psc ],
		              ax_hist_ylim,
		              color = dot_cols[ 0 ],
		              linestyle = "--"
		            )


#			ax.hold( True )

		plt.title( roi_name.upper() )

#		ylims.append( np.max( np.abs( ax.get_ylim() ) ) )

#		axes.append( ax )

#		ylim = np.max( ylims )

#		for ax in axes:

#			xlim = ax.get_xlim()

#			ax.plot( xlim, [ 0, 0 ], "k--" )

#			ax.set_ylim( ( -ylim, ylim  ) )
#			ax.set_xlim( xlim )


		plt.savefig( "/home/dmannion/im_temp/ns_ap_%s_%s.png" % ( conf[ "subj" ][ "subj_id" ], roi_name ) )
#	plt.show()


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
