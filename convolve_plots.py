import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ecco_v4_py as ecco
import numpy as np

basins = ['atl','pac','ind'] 
ecco_vars=['ADJsaltsurf', 'ADJthetasurf', 'adxx_qnet', 'adxx_tauu', 'adxx_tauv']
std_names=['eccoexf','cmip5','cmip6']
var_labels=dict(zip(ecco_vars,['SSS','SST','Q_{net}','\tau_u','\tau_v']))
basin_name=dict(zip(basins,['Atlantic Sector','Pacific Sector','Indian Sector']))
std_labels=['$ECCO_{std}$',r'$\langle ECCO_{std} \rangle$','$\overline{CMIP5_{std}}$','$CMIP6_{std}$']
std_labels_3d=['$ECCO_{std}$','$CMIP5_{std}$','$CMIP6_{std}$']
var_units = dict(zip(ecco_vars,['psu','$^{\circ}$C','W/m$^2$','N/m$^2$','N/m$^2$']))


def timeseries_plot(dataout_lib,cmip_ecco_map,ecco3d_ecco_map,plots_dir,cvar_ind=0):
    fig=plt.figure(figsize=(18,10)) 

    ylims={'atl':[0,2500],'pac':[0,1500],'ind':[0 ,2000]}

    pp=[[],[],[],[]]

    ip=1     
    for basin in basins:
        for var in ecco_vars:  
            
            dplot=[]
            for std_var in list(dataout_lib.keys()):
                if std_var=='eccoexf':
                    for cvar in ecco3d_ecco_map[var]:
                        dplot.append(dataout_lib[std_var]['{}X{}_absmean'.format(var,cvar)].mean(dim='year'))
                else:                    
                    cvar=cmip_ecco_map[var][cvar_ind]
                    dplot.append(dataout_lib[std_var]['{}X{}_absmean'.format(var,cvar)].mean(dim='year'))

            if ip>1:
                ax_old=ax
            ax=plt.subplot(3,5,ip)            
            for ipp in range(0,4):
                pp[ipp]=dplot[ipp].sel(basin=basin).plot(x='lag_years',linewidth=2,ax=ax,label=std_labels[ipp])

            if var=='ADJsaltsurf':
                ax.set_ylabel('deg C')
                ax.annotate( basin_name[basin],(0.1,0.8),xycoords='axes fraction',fontsize=14,fontweight='bold')
            else:
                ax.set_ylabel('')
            if basin=='ind':
                ax.set_xlabel('Lag (years)')
            else:
                ax.set_xlabel('')
            if basin=='atl':
                ax.set_title(var_labels[var],fontweight='bold')
            else:
                ax.set_title('')
            ax.set_ylim(ylims[basin])
            ip=ip+1
    ax.legend()        
    plt.suptitle('Impact of CMIP inter-model uncertainty vs in-model variability in surface variables on Southern Ocean MWFR heat content',fontweight='bold',y=0.94)
    plt.savefig(plots_dir+'/convolvetseries_CMIPECCO.png') 
    return fig

def plot_array_proj(subplotcode,x,y,dplot,cmin,cmax,cmap,lon0,lat0,lat_cut,paral,units,show_colorbar=True):
# Function to plot non-ecco data similarly to ecco-v4-py

    cbar=None
    ax = plt.subplot(subplotcode,projection=ccrs.LambertConformal(central_longitude=lon0,
                                        central_latitude=lat0,cutoff=lat_cut,standard_parallels=paral))
    out=ecco.plot_global(x,y,dplot,4326,cmin=cmin,cmax=cmax,ax=ax,\
                  plot_type='pcolormesh',show_colorbar=False,cmap=cmap,\
                  show_grid_lines = False)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(cmin,cmax))
    sm._A = []
    if show_colorbar:
        cbar_ax = plt.gcf().add_axes([0.94, 0.55-0.5*(subplotcode-233)/3, 0.02, 0.4])
        cbar = plt.colorbar(sm,ax=ax,cax=cbar_ax,label='[{}]'.format(units)) 
    ax.add_feature(cfeature.LAND, zorder=100)
    ax.add_feature(cfeature.COASTLINE,linewidth=0.5,zorder=101)
    ax.gridlines(crs=ccrs.PlateCarree(), 
                 linewidth=1.0,
                 color='black', 	
                 alpha=0.5, 
                 linestyle='--', 
                 draw_labels = False,zorder=102)
    return ax,cbar

def spatial_plot(var,basin,dataout_lib,ecco_exf_data,cmip5_spreads,cmip6_spreads,cmip_ecco_map,ecco3d_ecco_map,ecco_grid,plots_dir,clims,cvar_ind=0) :
       
    # Setup grid projection variables
    proj='LambertConformal'
    lon0={'atl':0,'pac':-160,'ind':80}
    lat0=-30
    lat_cut=-20
    paral=(-30,-50)

    dplot_1d=[]  
    dplot_3d=[]   
    for std_var in list(dataout_lib.keys()):
        if std_var=='eccoexf':   
            for ic,cvar in enumerate(ecco3d_ecco_map[var]):                
                dplot_1d.append(dataout_lib[std_var]['{}X{}_absmean'.format(var,cvar)].mean(dim='year'))     
                if ic==0:
                    dplot_3d.append(dataout_lib[std_var]['{}X{}'.format(var,cvar)].mean(dim='year').sel(basin=basin))                     
        else:
            cvar=cmip_ecco_map[var][cvar_ind]
            dplot_1d.append(dataout_lib[std_var]['{}X{}_absmean'.format(var,cvar)].mean(dim='year'))
            dplot_3d.append(dataout_lib[std_var]['{}X{}'.format(var,cvar)].mean(dim='year').sel(basin=basin))     

    # Find peak in basin integrated impacts
    dmax=0
    for ipp in range(0,4):
        if dplot_1d[ipp].sel(basin=basin).max()>dmax:
            imax=dplot_1d[ipp].sel(basin=basin).squeeze().argmax().data.compute()
            dmax=dplot_1d[ipp].sel(basin=basin).max()
            tmax=dplot_1d[ipp]['lag_years'][imax].data

    fig=plt.figure(figsize=(13,7))

    cmap='Reds'
    cmin,cmax=clims
    
    cvar=ecco3d_ecco_map[var][0]          
    dplot=ecco_exf_data[cvar]          
    f, ax=ecco.plot_proj_to_latlon_grid(ecco_grid.XC,ecco_grid.YC,dplot.compute(),user_lon_0=lon0[basin],user_lat_0=lat0,\
                                    projection_type=proj,plot_type='pcolormesh', show_colorbar=False,dx=1,dy=1,\
                                    cmap=cmap,subplot_grid=[2,3,1],cmin=cmin,cmax=cmax,lat_lim=lat_cut,
                                    parallels=paral)[:2]
    ax.set_title('${}$ {}'.format(var_labels[var],std_labels_3d[0]),fontweight='bold')

    cvar=cmip_ecco_map[var][0]    
    dplot=cmip5_spreads[cvar].mean(dim='time')
    ax,cbar=plot_array_proj(232,dplot.lon,dplot.lat,dplot,cmin,cmax,cmap,
                                lon0[basin],lat0,lat_cut,paral,var_units[var],show_colorbar=False)
    ax.set_title('${}$ {}'.format(var_labels[var],std_labels_3d[1]),fontweight='bold')

    dplot=cmip6_spreads[cvar].mean(dim='time')
    ax,cbar=plot_array_proj(233,dplot.lon,dplot.lat,dplot.T,cmin,cmax,cmap,
                                lon0[basin],lat0,lat_cut,paral,var_units[var])
    ax.set_title('${}$ {}'.format(var_labels[var],std_labels_3d[2]),fontweight='bold')

    cmax=20e3
    cmin=-cmax
    cmap='RdBu_r'      

    for idp,dplot in enumerate(dplot_3d):        
        dplot=dplot.where(np.isfinite(dplot.lag_years),drop=True).sel(lag_years=tmax,method='nearest').squeeze()

        if ('i' in dplot.dims) or ('i_g' in dplot.dims):
            f, ax=ecco.plot_proj_to_latlon_grid(ecco_grid.XC,ecco_grid.YC,dplot,user_lon_0=lon0[basin],user_lat_0=lat0,\
                                    projection_type=proj,plot_type='pcolormesh', show_colorbar=False,dx=1,dy=1,\
                                    cmap=cmap,subplot_grid=[2,3,4+idp],cmin=cmin,cmax=cmax,lat_lim=lat_cut,
                                    parallels=paral)[:2]
        else:
            if idp==1:
                scb=False
            else:
                scb=True
            ax,cbar=plot_array_proj(234+idp,dplot.lon,dplot.lat,dplot.transpose(),cmin,cmax,cmap,
                            lon0[basin],lat0,lat_cut,paral,'$^{\circ}$C',show_colorbar=scb)
        #print(ax.get_extent())
        ax.set_title('$\partial J/\partial {}$ * {}'.format(var_labels[var],std_labels_3d[idp]),fontweight='bold')

    plt.suptitle('${}$ variability in ECCO and between CMIP models,\n and mean projected impact on {} MWFR heat content from year {:2.1f} '.format(var_labels[var],basin_name[basin],tmax)
                 ,fontweight='bold')
    #plt.show()
    plt.savefig(plots_dir+var+'_'+basin+'_absall.png')