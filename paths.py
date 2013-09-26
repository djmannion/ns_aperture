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
	paths.exp_log = _get_exp_log_paths( conf, paths )
	paths.mvpa = _get_mvpa_paths( conf, paths )

	# add the filtered func type
	paths.func.filts = [ orig + orig.file().replace( "orig", "surf-filt" )
	                     for orig in paths.func.origs
	                   ]

	# add the mask
	paths.summ.mask = ( paths.summ.base +
	                    paths.summ.mean.file().replace( "mean", "mask" )
	                  )

	# add the group spec to the reg
	paths.reg.std_spec = fmri_tools.paths.Path( paths.reg.spec.dir(),
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

	loc.time_files = loc.base + ( file_base + "time_files" )
	loc.glm = loc.base + ( file_base + "glm" )
	loc.beta = loc.base + ( file_base + "beta" )
	loc.mask = loc.base + ( file_base + "mask" )

	return loc


def _get_exp_log_paths( conf, paths ):
	"""Get the paths for the logfiles"""

	log = fmri_tools.paths.PathsHandler()

	subj_id = conf.subj.subj_id
	exp_id = conf.exp.id

	file_base = "{subj_id:s}_{exp_id:s}".format( subj_id = subj_id, exp_id = exp_id )

	log.base = paths.base / "logs"

	log.seq = log.base + ( file_base + "_fmri_seq" )
	log.task = log.base + ( file_base + "_fmri_task" )

	return log


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

	ana.mask = ana.base + ( file_base + "mask" )

	ana.coef = ana.base + ( file_base + "coef" )
	ana.snr = ana.base + ( file_base + "snr" )

	ana.clust = ana.base + ( file_base + "clust" )

	return ana


def _get_mvpa_paths( conf, paths ):
	"""Get the paths for the MVPA analysis"""

	mvpa = fmri_tools.paths.PathsHandler()

	mvpa.base = paths.base / "mvpa"

	subj_id = conf.subj.subj_id
	exp_id = conf.exp.id

	file_base = "{subj_id:s}_{exp_id:s}-".format( subj_id = subj_id, exp_id = exp_id )

	mvpa.nodes = mvpa.base + ( file_base + "nodes" )
	mvpa.data = mvpa.base + ( file_base + "data" )
	mvpa.cond_info = mvpa.base + ( file_base + "cond_info" )
	mvpa.acc = mvpa.base + ( file_base + "acc" )

	return mvpa


def get_group_paths( conf ):

	grp = fmri_tools.paths.PathsHandler()

	grp.base = fmri_tools.paths.Path( "/labs/olmanlab/Data7T/NatSceneAperture/group_data" )

	file_base = "ns_aperture-"

	grp.log = grp.base + ( file_base + "log.log" )

	grp.avg = grp.base / "ns_aperture_avg"

	grp.coef = ( grp.base / "ret" ) + ( file_base + "coef" )
	grp.angle = ( grp.base / "ret" ) + ( file_base + "angle" )
	grp.vis_loc = ( grp.base / "ret" ) + ( file_base + "vis_loc" )

	grp.coh_test = grp.base + ( file_base + "coh_test" )
	grp.coh_clust = grp.base + ( file_base + "coh_clust" )

	grp.surf_mask = grp.base + ( file_base + "mask" )
	grp.rep_vol_mask = grp.base + ( file_base + "rep_vol_mask" )

	grp.clust_sim = grp.base / "cluster_sim" + "sim"
	grp.clust_script = grp.base / "cluster_sim" + "script"

	grp.loc = grp.base / "loc"
	grp.loc_test = grp.loc + ( file_base + "loc_test" )
	grp.loc_clust = grp.loc + ( file_base + "loc_clust" )

	grp.sl_info = ( grp.base / "mvpa" ) + ( file_base + "sl_info" )
	grp.sl_seed = ( grp.base / "mvpa" ) + ( file_base + "sl_seed" )
	grp.sl_disk = ( grp.base / "mvpa" ) + ( file_base + "sl_disk" )

	grp.acc = ( grp.base / "mvpa" ) + ( file_base + "acc" )
	grp.acc_clust = ( grp.base / "mvpa" ) + ( file_base + "acc_clust" )


	return grp




