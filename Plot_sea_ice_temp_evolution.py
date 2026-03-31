import xarray as xr
import matplotlib.pyplot as plt
import glob
import warnings
import tqdm
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import root_scalar
import pandas as pd
from matplotlib.legend_handler import HandlerTuple


#####

file_path = ''

#### Load data

warnings.filterwarnings('ignore')

def running_mean_years(T, window=5):
    w = np.ones(window) / window
    return np.apply_along_axis(
        lambda x: np.convolve(x, w, mode="same"),
        axis=0,
        arr=T
    )
colordict={}

temp_mean={}
T_smooth={}
temp_std={}

for scenario in ['ssp245', 'ssp370']:

    temp_sum = np.zeros((251,12,3))
    model_list_ssp = []
    model_list_hist = []
    for infile in tqdm.tqdm(sorted(glob.glob(file_path+'sitemptop/sitemptop*'+scenario+'*weighted*.nc'))):
        model=infile.rsplit('AO_',1)[1].split('_',1)[0]
        model_list_ssp.append(model)
    for infile in tqdm.tqdm(sorted(glob.glob(file_path+'sitemptop/sitemptop*historical*weighted*.nc'))):
        model=infile.rsplit('AO_',1)[1].split('_',1)[0]
        model_list_hist.append(model)
    
    for infile in tqdm.tqdm(sorted(glob.glob(file_path+'sitemptop/sitemptop*historical*weighted*.nc'))):
        model=infile.rsplit('AO_',1)[1].split('_',1)[0]
        if model in model_list_ssp and model in model_list_hist:
            ds=xr.load_dataset(infile)
            
            for n, year in enumerate(range(1850,2015)):
                for m, month in enumerate(range(1,13)):
                    try: # for some reasons, variable renaming did not always work, so sometimes siconc as variable name superseeded sitemptop, but contains temperature values
                        if np.isfinite(ds.sel(time=ds.time.dt.year.isin([year])).sitemptop[m]):
                            temp_sum[n,m,0] += ds.sel(time=ds.time.dt.year.isin([year])).sitemptop[m]
                            temp_sum[n,m,1] += 1
                            temp_sum[n,m,2] += (ds.sel(time=ds.time.dt.year.isin([year])).sitemptop[m])**2
                    except:
                        if np.isfinite(ds.sel(time=ds.time.dt.year.isin([year])).siconc[m]):
                            temp_sum[n,m,0] += ds.sel(time=ds.time.dt.year.isin([year])).siconc[m]
                            temp_sum[n,m,1] += 1
                            temp_sum[n,m,2] += (ds.sel(time=ds.time.dt.year.isin([year])).siconc[m])**2
                    
    for infile in tqdm.tqdm(sorted(glob.glob(file_path+'sitemptop/sitemptop*'+scenario+'*weighted*.nc'))):
        model=infile.rsplit('AO_',1)[1].split('_',1)[0]
        if model in model_list_ssp and model in model_list_hist:
            ds=xr.load_dataset(infile)
            for n, year in enumerate(range(2015,2100)):
                n += 165
                for m, month in enumerate(range(1,13)):
                    try: # for some reasons, variable renaming did not always work, so sometimes siconc as variable name superseeded sitemptop, but contains temperature values
                        if np.isfinite(ds.sel(time=ds.time.dt.year.isin([year])).sitemptop[m]):
                            temp_sum[n,m,0] += ds.sel(time=ds.time.dt.year.isin([year])).sitemptop[m]
                            temp_sum[n,m,1] += 1
                            temp_sum[n,m,2] += (ds.sel(time=ds.time.dt.year.isin([year])).sitemptop[m])**2
                    except:
                        if np.isfinite(ds.sel(time=ds.time.dt.year.isin([year])).siconc[m]):
                            temp_sum[n,m,0] += ds.sel(time=ds.time.dt.year.isin([year])).siconc[m]
                            temp_sum[n,m,1] += 1
                            temp_sum[n,m,2] += (ds.sel(time=ds.time.dt.year.isin([year])).siconc[m])**2
                        
    temp_mean[scenario]=(temp_sum[:,:,0] / temp_sum[:,:,1])-273.15
    temp_std[scenario]=np.sqrt((1/temp_sum[:,:,1])*(temp_sum[:,:,2]-(temp_sum[:,:,0]**2/temp_sum[:,:,1])))

    temp_mean[scenario] = np.hstack((temp_mean[scenario],temp_mean[scenario][:,[0]]))
    temp_std[scenario] = np.hstack((temp_std[scenario],temp_std[scenario][:,[0]]))
    
    T_smooth[scenario] = running_mean_years(temp_mean[scenario], window=7)
    
## load DMI satellite data
obs_label='SEAICE_ARC_PHY_CLIMATE_L4_MY_011_016'
DMI_df = pd.read_csv(file_path+'data/Final_DMI_IST_mean_monthlyfirst_and_SIC15.txt', index_col='year')

### plot


def fmt(x):
    s = f"{x:.1f}"
    return rf"{s} °C" if plt.rcParams["text.usetex"] else f"{s} °C"
fig1,ax1 = plt.subplots(1,1, figsize=(9,7), layout='constrained')
fig1.suptitle('', fontsize=17)
ice_bot=1.8
#ax1.invert_yaxis()
con = ax1.contourf((temp_mean['ssp370'][100:-1,:]-ice_bot).T/2., extent=(1950, 2100, 1,13), cmap='coolwarm', levels=np.linspace(-16,0,40,endpoint=True))
cons1 = ax1.contour((T_smooth['ssp370'][100:-1,:7]-ice_bot).T/2., extent=(1950, 2100, 1, 7), levels=sorted((T_smooth['ssp370'][100,[3,6]]-ice_bot)/2.),colors=['indigo','orange'], label='Model Mean')
cons2 = ax1.contour((T_smooth['ssp370'][100:-1,7:]-ice_bot).T/2., extent=(1950, 2100, 8, 13), levels=sorted((T_smooth['ssp370'][100,[9,11]]-ice_bot)/2.),colors=['indigo','firebrick'])
#cons3 = ax1.contour((T_smooth['ssp370'][100:-1,:]-ice_bot).T/2., extent=(1950, 2100, 1, 13), levels=sorted([np.mean((T_smooth['ssp370'][100,season_ind]-ice_bot)/2.) for season_ind in [
#    [11,0,1,2],[3,4,5],[6,7,8],[9,10]]]),colors=['black'])

fmt = {}
strs = ['winter\n-\nspring', 'spring\n-\nsummer', 'autumn\n-\nwinter', 'summer\n-\nautumn']
for l, s in zip(cons1.levels.tolist()+cons2.levels.tolist(), strs):
    fmt[l] = s
# Label every other level using strings
ax1.clabel(cons1, cons1.levels, manual=zip([1971,1962], [3,6]), fmt=fmt, fontsize=14)
ax1.clabel(cons2, cons2.levels, manual=zip([1971,1962], [9,11]), fmt=fmt, fontsize=14)
#ax1.clabel(cons1, [cons1.levels.tolist()[1]], manual=zip([2075], [6]), fmt=fmt, fontsize=14)
#ax1.clabel(cons2, [cons2.levels.tolist()[1]], manual=zip([2080], [2]), fmt=fmt, fontsize=14)

ax1.set_yticks(range(1,13))
ax1.set_yticklabels(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
ax1.set_xlabel('Years from 2026', fontsize=13)
ax1.set_xticks(range(1951,2077,25))
ax1.set_xticklabels([str(dy)+'y' for dy in range(-75,0,25)] + ['2026'] + ['+'+str(dy)+'y' for dy in range(25,51,25)])
ax1.invert_yaxis()
cons1 = ax1.contour((T_smooth['ssp245'][164:-1,:7]-ice_bot).T/2., extent=(2014, 2100, 1, 7), levels=sorted((T_smooth['ssp370'][100,[3,6]]-ice_bot)/2.),colors=['indigo','orange'], label='Model Mean', linestyles=':')
cons2 = ax1.contour((T_smooth['ssp245'][164:-1,7:]-ice_bot).T/2., extent=(2014, 2100, 8, 13), levels=sorted((T_smooth['ssp370'][100,[9,11]]-ice_bot)/2.),colors=['indigo','firebrick'], linestyles=':')


DMI_smooth = running_mean_years(DMI_df, window=5)
cons1 = ax1.contour((DMI_df.iloc[:,:9]-273.15-ice_bot).T/2., extent=(1982, 2023, 1, 9), levels=sorted((DMI_smooth[2,[3,6]]-273.15-ice_bot)/2.),colors=['white'], label='Observations', linestyles='-')
cons2 = ax1.contour((DMI_df.iloc[:,[8,9,10,11,0]]-273.15-ice_bot).T/2., extent=(1982, 2023, 9, 13), levels=sorted((DMI_smooth[2,[9,11]]-273.15-ice_bot)/2.),colors=['white'], linestyles='-')

line1, = ax1.plot([],[],ls='-',color='indigo',label='Model mean high emissions')
line2, = ax1.plot([],[],ls='-',color='orange')
line3, = ax1.plot([],[],ls='-',color='firebrick')
line4, = ax1.plot([],[],ls=':',color='indigo',label='Model mean intermediate emissions')
line5, = ax1.plot([],[],ls=':',color='orange')
line6, = ax1.plot([],[],ls=':',color='firebrick')
line7, = ax1.plot([],[],ls='-',color='white',label='Observations')
legend=ax1.legend([(line1, line2, line3), (line4, line5, line6), line7],
    ['Model mean\n(hist.+ssp370)\nref. 1950', 'Model mean\n(ssp245)\nref. 1950', 'Observation\n(DMI IST, relative)\nref. 1985'],
    handler_map={tuple: HandlerTuple(ndivide=None)}, title='Isotherms of seasonal transitions:', fancybox=True, ncol=3, framealpha=0.6, fontsize=12,
                 loc='upper center')
plt.setp(legend.get_title(),fontsize=13)
ax1.tick_params('both', labelsize=13, top=True, bottom=True,labeltop=False, labelbottom=True)
cbar = fig1.colorbar(con)
cbar.ax.set_ylabel('Sea ice temperature °C\nModel mean (CMIP6 multi-model mean hist+ssp370)',fontsize=13)
cbar.set_ticks(np.arange(-16,1,2))
cbar.ax.tick_params(labelsize=13)
ax1.set_title('Pan-Arctic Monthly Mean Sea Ice Temperature\n(SIC>15%)', fontsize=13)
fig1.savefig(file_path+'plots/IST_month_separation_Hovmoller_DMI_historical_ssps.png')




