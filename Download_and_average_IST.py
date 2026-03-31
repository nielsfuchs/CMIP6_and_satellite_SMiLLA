import copernicusmarine
import os
import numpy as np
import tqdm
import pandas as pd

# for Python script
os.environ['COPERNICUSMARINE_SERVICE_USERNAME'] = 'username'
os.environ['COPERNICUSMARINE_SERVICE_PASSWORD'] = 'fill_password_here'

copernicusmarine.login(username='nfuchs', password='fill_password_here')

# output array
mean = np.zeros((len(range(1982,2024)),13,2))

## load partial data and write to output array
for n_y, year in tqdm.tqdm(enumerate(range(1982,2024))):
    # Set parameters
    data_request = {
       "dataset_id_ist_dmi" : "cmems_obs_si_arc_phy_my_L4-DMIOI_P1D-m",
       "longitude" : [-179.97500610351562, 179.97500610351562], 
       "latitude" : [60, 89.94999694824219],
       "time" : [str(year)+"-01-01", str(year)+"-12-31"],
       "variables" : ["analysed_st", "mask", "sea_ice_fraction"]
    }

    # Load xarray dataset
    ds = copernicusmarine.open_dataset(
        dataset_id = data_request["dataset_id_ist_dmi"],
        minimum_latitude = data_request["latitude"][0],
        start_datetime = data_request["time"][0],
        end_datetime = data_request["time"][1],
        variables = data_request["variables"],
        coordinates_selection_method="strict-inside",
        dataset_version="202105"
    )

    if 'area_weights' not in locals():
        # Assuming your dataset has 'lat' and 'lon' as centers, 
        # you need to estimate bounds (e.g., halfway between centers)
        dlat = np.abs(ds.latitude[1] - ds.latitude[0])
        dlon = np.abs(ds.longitude[1] - ds.longitude[0])

        # Simple bound approximation
        lat_top = ds.latitude + dlat/2
        lat_bot = ds.latitude - dlat/2

        # Area in m^2 (using radian conversion)
        r_earth = 6371000
        area = (r_earth**2) * np.deg2rad(dlon) * (np.sin(np.deg2rad(lat_top)) - np.sin(np.deg2rad(lat_bot)))
        area_weights = np.repeat(abs(area.values)[:,np.newaxis],len(ds.longitude),axis=1) # Ensure positive area
        area_weights = area_weights/np.sum(area_weights)*len(ds.longitude)*len(ds.latitude)
        
    ds['area_weights'] = (('latitude','longitude'), area_weights)
    
    # SIC>15%
    mean[n_y,:-1,0] = (ds.analysed_st.where(ds.mask>6).where(ds.sea_ice_fraction>0.15).resample(time='1ME').mean().\
    weighted((ds['sea_ice_fraction'].where(ds.sea_ice_fraction>0.15).resample(time='1ME').mean().fillna(0))*ds['area_weights'])).mean(dim=['latitude','longitude']) # 
    
    mean[n_y,-1,0] = ds.analysed_st.where(ds.mask>6).where(ds.sea_ice_fraction>0.15).weighted((ds['sea_ice_fraction'].where(ds.sea_ice_fraction>0.15).fillna(0))*ds['area_weights']).mean(dim=['latitude','longitude','time'])
    
    # all
    mean[n_y,:-1,1] = (ds.analysed_st.where(ds.mask>6).resample(time='1ME').mean().\
    weighted((ds['sea_ice_fraction'].resample(time='1ME').mean().fillna(0))*ds['area_weights'])).mean(dim=['latitude','longitude'])
    
    mean[n_y,-1,1] = ds.analysed_st.where(ds.mask>6).weighted((ds['sea_ice_fraction'].fillna(0))*ds['area_weights']).mean(dim=['latitude','longitude','time'])
    del(ds)
    
    # some backup in between
    if n_y%5==0:
        np.save('backup_DMI_IST_mean_monthlyfirst_and_SIC15.npy', mean)

# output file with SIC>15 for plotting
df = pd.DataFrame(data=mean[:,:,0], index=range(1982,2024), columns=[str(i) for i in range(1,13)]+['yearmean'],copy=True)
df.to_csv('Final_DMI_IST_mean_monthlyfirst_and_SIC15.txt', index_label='year')
