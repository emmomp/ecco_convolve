# Convolving ECCO adjoint sensitivities with other fields
E Boland 2021 - contact emmomp@bas.ac.uk for any issues

Example notebooks that take ECCOv4 generated adjoint sensitivity fields (see [Forget et al. 2015](https://gmd.copernicus.org/articles/8/3071/2015/gmd-8-3071-2015.pdf)) and convolve with other fields.

Examples include convolving with ECCOv4 surface forcing standard deviation fields and inter-model standard deviations of surface fields from CMIP models (generated as in the [CMIP_ranges](https://github.com/emmomp/cmip_ranges) notebooks) 

The convolve_examples notebook loads ECCOv4 standard deviations, CMIP5 and CMIP6 intermodel standard deviations of various surface fields, convolves them with adjoint sensitivities[^1], and makes some example plots.

I recommend using the convolve_examples notebook as a starting point for generating and testing your own scripts that can be submitted to JASMIN's analysis machines.

Requires the ecco_v4_py library to read/plot ECCO data on the LLC grid, available here: https://github.com/ECCO-GROUP/ECCOv4-py/

[^1]: Held in the JASMIN ORCHESTRA workspace - see [here](https://help.jasmin.ac.uk/article/199-introduction-to-group-workspaces) for details of requesting access
