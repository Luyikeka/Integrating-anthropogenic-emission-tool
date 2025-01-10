import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import xesmf as xe
import pandas as pd

######=========providing the resolutions for input data and output data =========#######

ds_out = xe.util.grid_2d(70.125, 150, 0.25, 10.125, 60, 0.25) # longitude, and latitude range and resolution
ds_fine = xe.util.grid_2d(70.125, 150.125, 0.1, 10.125, 60.125, 0.1)
ds_coarse = xe.util.grid_2d(70.25, 150.25, 0.5, 10.25, 60.25, 0.5)
lon_50 = np.arange(70.25, 150.25, 0.5, dtype=np.float32)
lat_50 = np.arange(10.25, 60.25, 0.5, dtype=np.float32)
lon_25 = np.arange(70.125, 150, 0.25, dtype=np.float32)
lat_25 = np.arange(10.125, 60, 0.25, dtype=np.float32)
lon_10 = np.arange(70.10, 150.10, 0.1, dtype=np.float32)
lat_10 = np.arange(10.10, 60.10, 0.1, dtype=np.float32)

mapper = pd.read_csv('/home/yjzhang/meic2wrf/add_bg_emis/integrated_all/Integrated_mapper.csv')
mapper = mapper.set_index('MEIC')

######=========calculation of grid area =========#######
def ll_area(lat,res):  #input : lat : np.array((200, 320))
    Re=6371.392 ###km
    X=Re*np.cos(lat*(np.pi/180))*(np.pi/180)*res  #np.array((200, 320))
    Y=Re*(np.pi/180)*res  #float
    return X*Y

area_fine = ll_area(ds_fine.lat.values,0.1)
area_coarse = ll_area(ds_coarse.lat.values,0.5)


######===============read data by xarray ================#############
#####regridding waste emission
def regrid_aggregation(t, meic_spec):
    mapper_ceds = mapper.dropna()
    specs_ceds = mapper_ceds.index.values
    ceds_spec = mapper.loc[meic_spec, 'CEDS']
    if meic_spec in specs_ceds:###????
        ceds_spec = mapper.loc[meic_spec,'CEDS']
        dowst_path = '/mnt/beegfs/user/yjzhang/emis_data/CEDS/2017/CEDS_Glb_0.5x0.5_anthro_'+ ceds_spec+'__monthly_2017.nc'
        DS_dowst = xr.open_dataset(dowst_path)#.sel(lat=lat_50, lon=lon_50, method="nearest")
        re_dowst = DS_dowst['waste'][t].sel(lat=lat_50, lon=lon_50, method="nearest")#.isel(lat=lat_50, lon=lon_50, method="nearest")
        ds_coarse['dowst'] = re_dowst*0.001*area_coarse*2678400*1000000  #####convert unit here
        regridder_fcs = xe.Regridder(ds_coarse, ds_out, method, periodic=True,reuse_weights=True)
        re_wst = regridder_fcs(ds_coarse['dowst'].values)
        par = mapper.loc[meic_spec, 'partition']
        M = mapper.loc[meic_spec, 'weight']
        V = mapper.loc[meic_spec, 'if VOC']
        #####consider if it is VOC
        if V == 'Y':
            wst = re_wst*par/M
        else:
            wst = re_wst
    else:
         wst = np.zeros((200,320), dtype='float32')

    #####regridding shipping emission
    htap_spec = meic_spec
    doshp_path = '/mnt/beegfs/user/yjzhang/emis_data/HTAP/gridded/2017/edgar_HTAPv3_2017_'+ htap_spec +'.nc'
    DS_doshp = xr.open_dataset(doshp_path)#.sel(lat=lat_50, lon=lon_50, method="nearest")
    re_doshp = DS_doshp['HTAPv3_5_3_Domestic_shipping'][t].sel(lat=lat_10, lon=lon_10, method="nearest")#.isel(lat=lat_50, lon=lon_50, method="nearest")
    ds_fine['doshp'] = re_doshp*area_fine
    regridder_ffn = xe.Regridder(ds_fine, ds_out, method, periodic=True,reuse_weights=True)
    re_shp = regridder_ffn(ds_fine['doshp'].values)
    M = mapper.loc[meic_spec, 'weight']
    V = mapper.loc[meic_spec, 'if VOC']
    if V == 'Y':
        shp = re_shp/M
    else:
        shp = re_shp

    #####regridding aviation emission
    htap_spec = meic_spec
    doavi_path = '/mnt/beegfs/user/yjzhang/emis_data/HTAP/gridded/2017/edgar_HTAPv3_2017_'+ htap_spec +'.nc'
    DS_doavi = xr.open_dataset(doavi_path)#.sel(lat=lat_50, lon=lon_50, method="nearest")
    re_doavi = DS_doavi['HTAPv3_2_1_Domestic_Aviation'][t].sel(lat=lat_10, lon=lon_10, method="nearest")#.isel(lat=lat_50, lon=lon_50, method="nearest")
    ds_fine['doavi'] = re_doavi*area_fine
    regridder_ffn = xe.Regridder(ds_fine, ds_out, method, periodic=True,reuse_weights=True)
    re_avi = regridder_ffn(ds_fine['doavi'].values)
    if V == 'Y':
        avi = re_avi/M
    else:
        avi = re_avi

    ######===============put every regridded sectoral emission into one netcdf files================#############
    myds=xr.Dataset(
        {"waste": (( "lat","lon"), wst),
         "shipping": (( "lat","lon"), shp),
         "aviation": (( "lat","lon"), avi),},
        coords = {'lon': lon_25,
                  'lat': lat_25})

    #==============Output as NETCDF files================# 
    save_dir = '/home/yjzhang/meic2wrf/add_bg_emis/integrated_all'
    output = save_dir + '/' + 'regridded_aggregated_sectors' + '201701_' + meic_spec + '.nc'
    myds.to_netcdf(output,format="NETCDF3_CLASSIC")