69d68
<                  #grid_fname="CABLE_UNSW_GSWP3_gridinfo_0.5x0.5.nc",
71,72d69
<                  #mask_fname="gswp3_landmask_nomissing.nc",
<                  #mask_fname="SE_AUS_gswp3_landmask_nomissing.nc",
89,90c86,87
<         #self.grid_fname = "SE_aus_veg_types_AWAP_grid.nc"
<         self.grid_fname = "SE_aus_veg_types_AWAP_plus_LAI_fper_grid.nc"
---
>         #self.grid_fname = "/short/w35/mm3972/cable/src/CABLE-AUX/offline/gridinfo_mmy_MD_elev_orig_std_avg-sand_mask.nc"
>         self.grid_fname = "gridinfo_AWAP_EBF.nc"
95d91
<         #self.mask_fname = os.path.join("land_sea_mask/%s" % (mask_fname))
248,249c244
<     #cable_src = "../../src/trunk/trunk/"
<     cable_src = "../../src/trunk_DESICA_PFTs/trunk_DESICA_PFTs/"
---
>     cable_src = "../../src/trunk/trunk/"
251d245
<     #spinup_end_yr = 2000
252a247
>     #spinup_end_yr = 2000
263a259
>         #walltime = "4:00:00"
269c265
<         walltime = "1:00:00"
---
>         walltime = "6:30:00"
