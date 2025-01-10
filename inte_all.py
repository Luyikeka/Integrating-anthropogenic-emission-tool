import pandas as pd
import rioxarray
import xarray as xr
import numpy as np
import geopandas as gpd
import fnmatch
import xarray as xr
import os
import matplotlib.pyplot as plt
import os

def ll_area(lat,res):  #input : lat : np.array((200, 320))
    Re=6371.392 ###km
    X=Re*np.cos(lat*(np.pi/180))*(np.pi/180)*res  #np.array((200, 320))
    Y=Re*(np.pi/180)*res  #float
    return X*Y

#def emis_union(ceds_dir,meic_dir,save_dir, spec_ceds, spec_meic):
def emis_union(ceds_dir,meic_dir,save_dir, spec_ceds, spec_meic,mon,mon_id,year):
    #####-------------------------------CAMs----------------------------------------####
    pre_ceds = 'CEDS_Glb_0.5x0.5_anthro_'
    post_ceds = '__monthly_'+year+'.nc'
    CEDS = ceds_dir + '/' + pre_ceds + spec_ceds + post_ceds
    if spec_ceds == 'BC': CEDS = ceds_dir + '/' + 'CEDS_Glb_0.5x0.5_anthro_' + 'BC__monthly_'+'2016'+'.nc'
    ds = rioxarray.open_rasterio(CEDS, masked=True)###read by rio, Easy to clip
    re_lat = np.arange(ds.y.values.min(),ds.y.values.max(), 0.25)###re interplote grids
    re_lon = np.arange(ds.x.values.min(),ds.x.values.max(), 0.25)
    emis_all = ds.interp(y=re_lat, x=re_lon)###dataset be interploted from 0.1 to 0.25
    emis_all = emis_all.isel(time=mon_id)
    lon_arange = np.arange(70.125, 150, 0.25, dtype=np.float32)##MEIC lat-lon ranges
    lat_arange = np.arange(10.125, 60, 0.25, dtype=np.float32)
    res_tomeic = emis_all.sel(x=lon_arange, y=lat_arange, method="nearest")####intercept meic square ###排放总量

    
    #########------------------------calculate area of meic grids 0.25 degree---------------########
    lon_2d, lat_2d = np.meshgrid(lon_arange, lat_arange) 
    area = ll_area(lat_2d,0.25)

    #########---------------convert units same as meic---------------########
    ####区分VOC species
    mapper = pd.read_csv('/home/yjzhang/meic2wrf/add_bg_emis/integrated_all/Integrated_mapper.csv')
    mapper = mapper.set_index('MEIC')
    meic_spec = spec_meic
    par = mapper.loc[meic_spec, 'partition']
    M = mapper.loc[meic_spec, 'weight']  
    V = mapper.loc[meic_spec, 'if VOC'] 
    if V == 'Y':
        unit_tomeic = res_tomeic*0.001*area*2678400*1000000*par/M
    else:
        unit_tomeic = res_tomeic*0.001*area*2678400*1000000*par

    #########---------------clip by shape file---------------########
    country = gpd.read_file(r'/home/yjzhang/shpfiles/World_GIS_data/country.shp')
    China_shp = country[country['CNTRY_NAME'] == 'China']#####台湾， 香港，澳门
    province = gpd.read_file(r'/home/yjzhang/shpfiles/全国数据shp/全国数据/分省.shp')
    Taiwan_shp = province[province['行政区划_c'] == '台湾省']
    ##Difference of Polygons
    mChina = gpd.overlay(China_shp, Taiwan_shp, how = 'difference')
    unit_tomeic.rio.write_crs("epsg:4326", inplace=True)
    ceds_clipped = unit_tomeic.rio.clip(mChina.geometry, mChina.crs, drop=False, invert=True)
    
    #####execute the regridding programmes
    #regrid_aggregation(t, meic_spec)
    
    ###read the ready intermediate aggregation emission data
    agg_path = '/mnt/beegfs/user/yjzhang/emission/integrated_all/agg_sectors/regridded_aggregated_sectors1'+year +mon+'_'+meic_spec+'.nc'
    ds_agg = rioxarray.open_rasterio(agg_path, masked=True)
    ds_agg.rio.write_crs("epsg:4326", inplace=True)
    agg_clipped = ds_agg.rio.clip(mChina.geometry, mChina.crs, drop=False, invert=True)
    DS_agg = xr.open_dataset(agg_path)

    
    ####masking waste emission
    ####读取CEDS数据需要考虑mapper table 中的空值
    dowst = DS_agg['waste'].values
    dowst_clip = np.nan_to_num(agg_clipped['waste'].values, nan = 0)
    dms_waste = dowst - dowst_clip[0][::-1]
    
    ####读取CEDS数据需要考虑mapper table 中的空值
    doagr = DS_agg['agriculture'].values
    doagr_clip = np.nan_to_num(agg_clipped['agriculture'].values, nan = 0)
    dms_agr = doagr - doagr_clip[0][::-1]

    ####masking shipping emission
    doshp = DS_agg['shipping'].values
    doshp_clip = np.nan_to_num(agg_clipped['shipping'].values, nan = 0)
    dms_shipping = doshp - doshp_clip[0][::-1]

    ####masking aviation emission
    all_avi = DS_agg['aviation'].values
    #doavi_clip = np.nan_to_num(agg_clipped['aviation'].values, nan = 0)
    #dms_aviation = doavi - doavi_clip[0][::-1]
    
    

    #########---------------read meic files---------------########
    ent_dir = meic_dir
    i = '*' + '_' + spec_meic+ '.' + '*'
    fn_act = ent_dir+'/' + fnmatch.filter(fnmatch.filter(os.listdir(ent_dir), i), '*agr*nc')[0]###the intermediate files mentioned above, include five sections as different variables in nc files
    fn_idt = ent_dir+'/' + fnmatch.filter(fnmatch.filter(os.listdir(ent_dir), i), '*ind*nc')[0]
    fn_pwr = ent_dir+'/' + fnmatch.filter(fnmatch.filter(os.listdir(ent_dir), i), '*pow*nc')[0]
    fn_rdt = ent_dir+'/' + fnmatch.filter(fnmatch.filter(os.listdir(ent_dir), i), '*res*nc')[0]
    fn_tpt = ent_dir+'/' + fnmatch.filter(fnmatch.filter(os.listdir(ent_dir), i), '*tra*nc')[0]
    f_act = xr.open_dataset(fn_act)
    f_idt = xr.open_dataset(fn_idt)
    f_pwr = xr.open_dataset(fn_pwr)
    f_rdt = xr.open_dataset(fn_rdt)
    f_tpt = xr.open_dataset(fn_tpt)
    act = f_act['z'][:].values.reshape((200, 320),)[::-1]
    act = np.where(act > 0.0, act*1, 0.0)
    idt = f_idt['z'][:].values.reshape((200, 320),)[::-1]
    idt = np.where(idt > 0.0, idt*1, 0.0)
    pwr = f_pwr['z'][:].values.reshape((200, 320),)[::-1]
    pwr = np.where(pwr > 0.0, pwr*1, 0.0)
    rdt = f_rdt['z'][:].values.reshape((200, 320),)[::-1]
    rdt = np.where(rdt > 0.0, rdt*1, 0.0)
    tpt = f_tpt['z'][:].values.reshape((200, 320),)[::-1]
    tpt = np.where(tpt > 0.0, tpt*1, 0.0)####用零取代缺失值 -9999
    
    ###---------------re calculate sectors---------------####Modify here to scaling sectors of CEDS to MEIC

    pwr_union = np.nan_to_num(ceds_clipped['energy'], nan = 0) + pwr
    res_union = np.nan_to_num(ceds_clipped['residential'], nan = 0) + np.nan_to_num(ceds_clipped['solvents'], nan = 0) + rdt
    idt_union = np.nan_to_num(ceds_clipped['industrial'], nan = 0) + idt
    shp_out = np.nan_to_num(ceds_clipped['ships'], nan = 0)
    shp_union = shp_out + dms_shipping
    swd_out = np.nan_to_num(ceds_clipped['waste'], nan = 0)#### to add domestic shipping emission 
    swd_union = swd_out + dms_waste
    #avi_out = np.nan_to_num(ceds_clipped['aviation'], nan = 0)
    avi_union = all_avi####aviation outside of China
    tpt_union = np.nan_to_num(ceds_clipped['transportation'], nan = 0) + tpt
    act_union = np.nan_to_num(ceds_clipped['agriculture'], nan = 0) + dms_agr
    sum_union = pwr_union + res_union + idt_union  + shp_union + swd_union + tpt_union + act_union

    #==============Create a xarray file to serve reunion emission data================# 
    myds=xr.Dataset(
        {"energy": (( "lat","lon"), pwr_union),
         "residential": (( "lat","lon"), res_union),
         "industry": (( "lat","lon"), idt_union),
         "agriculture": (( "lat","lon"), act_union),
         "transportation": (( "lat","lon"), tpt_union),
         "waste": (( "lat","lon"), swd_union),
         "shipping": (( "lat","lon"), shp_union),
         "aviation": (( "lat","lon"), avi_union),
         "sum": (( "lat","lon"), sum_union),},
        coords = {'lon': lon_arange,
                  'lat': lat_arange})
    if V == 'Y':
        myds.attrs['unit'] = 'ton/month/grid'
    else:
        myds.attrs['unit'] = 'million ton/month/grid'

    myds.attrs['conventions'] = 'NETCDF3_CLASSIC'
    myds.attrs['comments'] = 'Integrated inventories include more sectoral emission from both global inventories and regional inventories and uniform VOC speciation as MOZART chemistry mechanism.'
    myds.attrs['projection'] = 'Latitude-Longitude gridded data at a 0.25 x 0.25 decimal degrees spatial resolution.'
    myds.attrs['authors'] = 'Files prepared by Yijuan Zhang: University of Bremen.'
    myds.attrs['source'] = 'Data are available at the '
    myds.attrs['title'] = 'Integrated anthropogenic emission inventory for China in ' + year
    #==============Output as NETCDF files================# 
    #output = '/mnt/beegfs/user/yjzhang/emission/reunion/' + 'union_' + '201707_' + i.split('*')[1] + '.nc' 
    output_spec = mapper.loc[meic_spec,'output species']
    output = save_dir + '/' + 'Integrated_Anthropogenic_CEDS-bg_' + year +'_'+mon+'_' + output_spec + '_0.25x0.25' +'.nc'#####automatically the date
    myds.to_netcdf(output,format="NETCDF3_CLASSIC")
    

    
####Iteration of emission data by month
ceds_dir = '/mnt/beegfs/user/yjzhang/emis_data/CEDS'
#meic_dir = '/mnt/beegfs/user/yjzhang/emission/meic/201701/merge_mozart/' ###check the path and data time before run code
save_dir = '/mnt/beegfs/user/yjzhang/emission/integrated_all/'

year = '2018'
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
months_index = [0,1,2,3,4,5,6,7,8,9,10,11]
#months = ['07','08','09','10','11','12']
#months_index = [6,7,8,9,10,11]
#months = ['01']
#months_index = [0]

mapper = pd.read_csv('/home/yjzhang/meic2wrf/add_bg_emis/integrated_all/Integrated_mapper.csv')
mapper = mapper.set_index('MEIC')
for mon, mon_id in zip(months, months_index):
    #meic_dir = '/mnt/beegfs/user/yjzhang/emission/meic/'+year+mon+'/merge_mozart/'
    meic_dir = '/mnt/beegfs/user/yjzhang/emission/meic/'+year+'/'
    save_dir = '/mnt/beegfs/user/yjzhang/emission/integrated_all/'+year + mon +'/'
    t = mon_id
    if not os.path.exists(save_dir): os.makedirs(save_dir)
    for meic_spec in mapper.index.values:
        spec_meic = meic_spec
        spec_ceds = mapper.loc[meic_spec, 'CEDS']
        emis_union(ceds_dir,meic_dir,save_dir, spec_ceds, spec_meic,mon,mon_id,year)
        print(meic_spec)