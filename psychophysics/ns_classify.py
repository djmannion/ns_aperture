import numpy as np

import psychopy, psychopy.visual, psychopy.core, psychopy.event
import psychopy.filters

import stimuli.utils

im_path = "../im_db"

sel_path = "../data/im_sel.npy"
sel_data = np.load( sel_path )

i_img = np.where( sel_data == 1 )[ 0 ]

ap_diam_deg = 4.0
ap_ecc_deg = 3.0

ap_diam_pix = ap_diam_deg * 60.0
ap_ecc_pix = ap_ecc_deg * 60.0

reg_rows = np.arange( ( 1024 / 2 ) - np.round( ap_diam_pix / 2 ),
                      ( 1024 / 2 ) + np.round( ap_diam_pix / 2 )
                    )
reg_cols_l = np.arange( ( 1536 / 2 ) - ap_ecc_pix - np.round( ap_diam_pix / 2 ),
                        ( 1536 / 2 ) - ap_ecc_pix + np.round( ap_diam_pix / 2 )
                      )

reg_cols_r = np.arange( ( 1536 / 2 ) + ap_ecc_pix - np.round( ap_diam_pix / 2 ),
                        ( 1536 / 2 ) + ap_ecc_pix + np.round( ap_diam_pix / 2 )
                      )

win = psychopy.visual.Window( ( 1536, 1024 ),
                              fullscr = False,
                              allowGUI = True
                            )

status = psychopy.visual.TextStim( win = win,
                                   text = "",
                                   pos = ( 0, -0.75 )
                                 )

mask = psychopy.filters.makeMask( 256,
                                  radius = ap_diam_pix / 256,
                                  shape = "circle",
                                  fringeWidth = 0.0
                                )

fixation = psychopy.visual.Circle( win = win, radius = 3, units = "pix", fillColor = ( -1, -1, -1, 1 ) )

max_cd = 60.0

keep_going = True

while keep_going:

	# pick a random image
	im_l_id = i_img[ np.random.randint( i_img.shape[ 0 ] ) ]

	img_l = stimuli.utils.read_van_hateren( im_l_id + 1,
	                                        path = im_path,
	                                        scale = "norm",
	                                        region = ( reg_rows, reg_cols_l ),
	                                        pad = True,
	                                        flip = True
	                                      )

	same = bool( np.random.randint( 2 ) )

#	same = True

	if not same:
		im_r_id = i_img[ np.random.randint( i_img.shape[ 0 ] ) ]
	else:
		im_r_id = im_l_id

	img_r = stimuli.utils.read_van_hateren( im_r_id + 1,
	                                        path = im_path,
	                                        scale = "norm",
	                                        region = ( reg_rows, reg_cols_r ),
	                                        pad = True,
	                                        flip = True
	                                      )

	imgs = [ psychopy.visual.PatchStim( win = win,
	                                    tex = im_tex,
	                                    size = im_tex.shape[ 0 ],
	                                    units = "pix",
	                                    mask = mask,
	                                    pos = ( im_pos, 0 )
	                                  )
	         for ( im_tex, im_pos ) in zip( ( img_l, img_r ), ( -ap_ecc_pix, +ap_ecc_pix ) )
	       ]

	any( im.draw() for im in imgs )

	status.setText( "Image %d, Image %d" % ( im_l_id + 1, im_r_id + 1 ) )

	status.draw()

	fixation.draw()

	win.flip()

	valid_key = False

	while not valid_key:

		# wait for a key response
		keys = psychopy.event.waitKeys()

		for key in keys:

			if key in [ "a", "l", "q" ]:
				valid_key = True

			if key == "q":
				keep_going = False


psychopy.core.quit()
