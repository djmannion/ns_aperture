import os, os.path

import fmri_tools.paths


def get_exp_paths( conf ):

	try:
		code_dir = os.environ[ "CODE_DIR" ]
	except KeyError:
		print "CODE_DIR not set in shell environment"
		raise

	exp = fmri_tools.paths.PathsHandler()

	exp.base_dir = os.path.join( code_dir, "ns_apertures" )

	return exp


def get_subj_paths( conf ):
	"""Get the path structure for a given subject"""

	base_dir = os.path.join( "/labs/olmanlab/Data7T/NatSceneAperture/subj_data",
	                         conf.subj.subj_id
	                       )

	paths = fmri_tools.paths.get_default_paths( base_dir = base_dir,
	                                            subj_id = conf.subj.subj_id,
	                                            study_id = conf.exp.id,
	                                            n_runs = conf.subj.n_runs
	                                          )

	paths.ana = _get_ana_paths( conf, paths )
	paths.loc = _get_loc_paths( conf, paths )

	# add the group spec to the reg
	paths.reg.group_spec = fmri_tools.paths.Path( paths.reg.spec.dir(),
	                                              "std.141." + conf.subj.subj_id
	                                            )

	return paths


def _get_loc_paths( conf, paths ):
	"""Get the paths for the localisers"""

	loc = fmri_tools.paths.PathsHandler()

	loc.base = paths.base / "loc_analysis"

	subj_id = conf.subj.subj_id
	exp_id = conf.exp.id + "_loc"

	file_base = "{subj_id:s}_{exp_id:s}-".format( subj_id = subj_id,
	                                              exp_id = exp_id
	                                            )

	loc.glm = loc.base + ( file_base + "glm" )
	loc.beta = loc.base + ( file_base + "beta" )
	loc.mask = loc.base + ( file_base + "mask" )

	return loc


def _get_ana_paths( conf, paths ):
	"""Get the paths for the analysis"""

	ana = fmri_tools.paths.PathsHandler()

	ana.base = paths.base / "analysis"

	subj_id = conf.subj.subj_id
	exp_id = conf.exp.id

	file_base = "{subj_id:s}_{exp_id:s}-".format( subj_id = subj_id, exp_id = exp_id )

	ana.stim_times = ana.base + ( file_base + "stim_times" )

	ana.glm = ana.base + ( file_base + "glm" )
	ana.beta = ana.base + ( file_base + "beta" )

	ana.bltc = ana.base + ( file_base + "bltc" )
	ana.bl = ana.base + ( file_base + "bl" )
	ana.psc = ana.base + ( file_base + "psc" )

	ana.vl = ana.base + ( file_base + "vis_loc_rois" )
	ana.mask = ana.base + ( file_base + "mask" )

	return ana



def get_group_paths( conf ):

	grp = fmri_tools.paths.PathsHandler()

	grp.base = fmri_tools.paths.Path( "/labs/olmanlab/Data7T/NatSceneApertures/group_data" )

	file_base = "ns_apertures"

	return grp




