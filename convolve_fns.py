#!/usr/bin/env python
# coding: utf-8

import numpy as np
import xarray as xr
import pyresample

def convolve_ecco(ecco_data,convolve_data,ecco_convolve_map,ecco_grid=None,regrid=True,dir_out=None):
    
    # Function to convolve ECCO adjoint sensitivities with given fields to produce convolutions
    # ECCO data regridded to 1deg grid automatically, and, if convolve data is time-dependent, the sensetivities are resampled in time to match
    # 
    # Parameters
    # ----------
    # ecco_data : xarray dataset 
    #      ECCO sensitivities to convolve
    # convolve_data : xarray dataset
    #      Fields to convolve with ECCO data
    # ecco_convolve_map : dict
    #      Mapping from each ecco_data variable to convolve_data variable or list of variables. Used to match variables to convolve together
    # regrid : logical
    #      If True (default), regrid ecco_data to convolve_data grid, must have variables 'lon' and 'lat' 
    # dir_out : str
    #      Directory to write data_out to. If None (default), no data written
    #
    # Returns
    # -------
    # data_all : xarray dataset
    #      Convolved data, returned as variables with names '[ecco_var]X[convolve_var]', 
    #      where 'ecco_var' and 'convolve_var' are variables from ecco_data and convolve_data respectively
    
    if regrid:   
        regridder=setup_regrid(ecco_grid.XC,ecco_grid.YC,convolve_data.lon,convolve_data.lat) 
     
    data_all = []
    
    for var in list(ecco_convolve_map.keys()):
        data_all_var=[]
        for basin in list(ecco_data.basin.values):
            data_all_basin=[]
            for year in list(ecco_data.year.values):
                print(var,basin,year)
                # Select sensitivity for basin and year
                sens = ecco_data[var].sel(basin=basin,year=year)
                # Regrid in space to match convolve_data if necessary
                if regrid: 
                    sens=repeat_regrid(sens,regridder,convolve_data.lon,convolve_data.lat)
                    if 'year' in ecco_data['time'].dims:
                        sens['time']=ecco_data['time'].sel(year=year)
                    else:
                        sens['time']=ecco_data['time']
                    
                data_out=None
                for cvar in ecco_convolve_map[var]: 
                    # Resample in time to match convolve_data if necessary
                    if ('time' in convolve_data) :
                        if convolve_data[cvar]['time'].size>1:
                            if 'time' not in sens.dims:
                                sens=sens.swap_dims({'lag_years':'time'})
                            sens = sens.interp(time=convolve_data[cvar]['time'])
                            data_convolve = sens*convolve_data[cvar] 
                            data_convolve=data_convolve.swap_dims({'time':'lag_years'})
                            data_convolve=data_convolve.dropna('lag_years',how='all')
                        else:
                            data_convolve = sens*convolve_data[cvar] 
                    else:
                    # Perform convolution 
                        data_convolve = sens*convolve_data[cvar] 
                                                
                    # Take sums (sensitivities already area weighted)
                    dims = list(data_convolve.dims)
                    dims.remove('lag_years')
                    data_convolve_sum=data_convolve.sum(dim=dims)
                    data_convolve_abssum=np.abs(data_convolve).sum(dim=dims)
                    
                    # Put in a dataset
                    out_dict={'{}X{}'.format(var,cvar):data_convolve,
                            '{}X{}_sum'.format(var,cvar):data_convolve_sum,
                            '{}X{}_abssum'.format(var,cvar):data_convolve_abssum}
                    if data_out is None:
                        data_out=xr.Dataset(data_vars=out_dict)
                    else:
                        data_out=data_out.assign(out_dict)
                        
                data_out['basin']=(('basin',[basin,]))
                data_out['year']=(('year',[year,]))   
                    
                data_all_basin.append(data_out)
            data_all_var.append(xr.concat(data_all_basin,'year'))
         
        data_all_var=xr.combine_nested(data_all_var,['basin'])
        data_all.append(data_all_var)
        # Write to file if required
        if dir_out is not None:
            data_all_var.to_netcdf(dir_out+'/ECCOconv_{}.nc'.format(var))  
            print('Written to',dir_out+'/ECCOconv_{}.nc'.format(var))
                
    data_all = xr.merge(data_all)
    return data_all

def setup_regrid(xc,yc,new_lon=np.linspace(-179,180,360),new_lat=np.linspace(-89.5,89.5,180)):
    orig_grid = pyresample.geometry.SwathDefinition(lons=xc.values.ravel(), lats=yc.values.ravel())
    yi,xi=np.meshgrid(new_lat,new_lon)
    new_grid  = pyresample.geometry.GridDefinition(lons=xi,lats=yi)
    resample_data = pyresample.kd_tree.get_neighbour_info(orig_grid, new_grid, 100000, neighbours=1)
    return resample_data   

def repeat_regrid(ds,resample_data,new_lon=np.linspace(-179,180,360), new_lat=np.linspace(-89.5,89.5,180),loop_dim='lag_years'):    
    grid_shape=[new_lon.size,new_lat.size]
    stack_dims=ds.dims[1:]
    foo = pyresample.kd_tree.get_sample_from_neighbour_info('nn', grid_shape, ds.stack(z=stack_dims).transpose(...,loop_dim).values,
                                              resample_data[0], resample_data[1],resample_data[2])    
    ds2=xr.DataArray(foo,dims=['lon','lat',loop_dim],coords={'lon':(('lon'),new_lon.data),'lat':(('lat'),new_lat.data),loop_dim:(loop_dim,ds[loop_dim].data)})
    return ds2


