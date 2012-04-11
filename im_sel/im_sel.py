import numpy as np

import psychopy, psychopy.visual, psychopy.core, psychopy.event
import psychopy.filters

import stimuli.utils

im_path = "../im_db"

im_details = stimuli.utils.read_van_hateren_settings( im_path )

n_images = im_details.shape[ 0 ]

sel_path = "../data/im_sel.npy"

try:
	sel_data = np.load( sel_path )
except IOError:
	sel_data = np.empty( n_images )
	sel_data.fill( np.NAN )

n_good_images = 150

# find the indices of the images still todo
im_todo = np.where( np.isnan( sel_data ) )[ 0 ]

im_mode = "classify"

# if there are no NaNs left
if im_todo.shape[ 0 ] == 0:

	print "All classifications complete. Switching to cull mode."

	im_mode = "cull"

	# get the number of 'ok' images
	im_todo = np.where( sel_data == 1 )[ 0 ]


# randomise (in-place)
np.random.shuffle( im_todo )

ap_diam_deg = 3.0
ap_ecc_deg = 3.0

ap_diam_pix = ap_diam_deg * 60.0
ap_ecc_pix = ap_ecc_deg * 60.0


win = psychopy.visual.Window( ( 1536, 1024 ),
                              fullscr = False,
                              allowGUI = True
                            )

aps = [ psychopy.visual.Circle( win = win,
                                radius = ap_diam_pix / 2.0,
                                units = "pix",
                                lineColor = ( 1.0, 0.0, 0.0 ),
                                pos = ( ap_ecc_pix * merid_dir, 0 )
                              )
       for merid_dir in ( -1, +1 )
      ]

img = stimuli.utils.read_van_hateren( 1,
                                      path = im_path,
                                      scale = "norm",
                                      pad = True,
                                      flip = True
                                    )

img_tex = psychopy.visual.PatchStim( win = win,
                                     tex = img,
                                     size = img.shape[ 0 ],
                                     units = "pix"
                                   )

max_cd = 60.0

keep_going = True

while keep_going:

	for i_im in im_todo:

		if not keep_going:
			break

		if np.sum( sel_data == 1 ) == n_good_images and im_mode == "cull":

			print "Good images at the right number. Exiting"

			keep_going = False
			break

		try:
			# load the image
			img = stimuli.utils.read_van_hateren( i_im + 1,
			                                      path = im_path,
			                                      scale = "norm",
			                                      pad = True,
			                                      flip = True
			                                    )
		except IOError:
			sel_data[ i_im ] = 0
			continue

		# update the texture
		img_tex.setTex( img )

		img_tex.draw()

		for ap in aps:
			ap.draw()

		win.flip()

		valid_key = False

		while not valid_key:

			# wait for a key response
			keys = psychopy.event.waitKeys()

			for key in keys:

				if key in [ "a", "l", "q" ]:
					valid_key = True

				if key == "a":
					sel_data[ i_im ] = 0
				elif key == "l":
					sel_data[ i_im ] = 1
				elif key == "q":
					keep_going = False

np.save( sel_path, sel_data )

psychopy.core.quit()
