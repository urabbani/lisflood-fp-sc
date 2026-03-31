
# coding: utf-8


import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
from matplotlib import colors

import sys
import logging


def gdal_readmap(file_name, file_format, give_geotrans=False):
    """ Read geographical file into memory
    Dependencies are osgeo.gdal and numpy
    Input:
        file_name: -- string: reference path to GDAL-compatible file
        file_format: -- string: file format according to GDAL acronym
        (see http://www.gdal.org/formats_list.html)
        give_geotrans (default=False): -- return the geotrans and amount of 
            cols/rows instead of x, y axis
    Output (if give_geotrans=False):
        x: -- 1D np-array: x-axis
        y: -- 1D np-array: y-axis
        data:           -- 2D np-array: raster data
        fill_val         -- float:       fill value
    Output (if give_geotrans=True):
        geotrans: -- 6-digit list with GDAL geotrans vector
        size: -- 2-digit tuple with (cols, rows)
        data:           -- 2D np-array: raster data
        fill_val         -- float:       fill value
    """
    # Open file for binary-reading
    mapFormat = gdal.GetDriverByName(file_format)
    mapFormat.Register()
    ds = gdal.Open(file_name)
    if ds is None:
        logging.warning('Could not open {:s} Shutting down'.format(file_name))
        sys.exit(1)
        # Retrieve geoTransform info
    geotrans = ds.GetGeoTransform()
    originX = geotrans[0]
    originY = geotrans[3]
    resX = geotrans[1]
    resY = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = np.linspace(originX+resX/2, originX+resX/2+resX*(cols-1), cols)
    y = np.linspace(originY+resY/2, originY+resY/2+resY*(rows-1), rows)
    # Retrieve raster
    RasterBand = ds.GetRasterBand(1)   # there's only 1 band, starting from 1
    data = RasterBand.ReadAsArray(0, 0, cols, rows)
    fill_val = RasterBand.GetNoDataValue()
    RasterBand = None
    ds = None
    if give_geotrans==True:
        return geotrans, (ds.RasterXSize, ds.RasterYSize), data, fill_val
        
    else:
        return x, y, data, fill_val



def read_axes(fn):
    """
    Retrieve x and y ax from a raster datasets (GDAL compatible)
    """
    ds = gdal.Open(fn)
    geotrans = ds.GetGeoTransform()
    originX = geotrans[0]
    originY = geotrans[3]
    resX = geotrans[1]
    resY = geotrans[5]
    cols = ds.RasterXSize
    rows = ds.RasterYSize
    x = np.linspace(originX+resX/2, originX+resX/2+resX*(cols-1), cols)
    y = np.linspace(originY+resY/2, originY+resY/2+resY*(rows-1), rows)
    ds = None
    return x, y


def contingency_map(array1, array2, threshold1=0., threshold2=0.):
    """
    Establish the contingency between array1 and array2.
    Returns an array where 
    1 means only array2 gives a value > threshold1, 
    2 means only array1 gives a values > threshold2,
    3 means array1 gives a value > threshold1, and array2 a value > threshold2
    0 means both arrays do not give a value > threshold1, 2 respectively
    
    function returns the threshold exceedance (0-1) of array 1 and 2, as well as the contingency map
    """
    array1_thres = array1 > threshold1
    array2_thres = array2 > threshold2
    contingency = array1*0
    contingency += np.int16(array2_thres)
    contingency += np.int16(array1_thres)*2
    return array1_thres, array2_thres, contingency


def hit_rate(array1, array2):
    """
    calculate the hit rate based upon 2 boolean maps. (i.e. where are both 1)
    """
    # count the number of cells that are flooded in both array1 and 2
    idx_both = np.sum(np.logical_and(array1, array2))
    idx_1 = np.sum(array1)
    return float(idx_both)/float(idx_1)


def false_alarm_rate(array1, array2):
    """
    calculate the false alarm rate based upon 2 boolean maps. (i.e. amount of cells where array2 is True but array1 False)
    """
    # count the number of cells that are flooded in both array1 and 2
    idx_2_only = np.sum(np.logical_and(array2, array1!=1))
    idx_2_total = np.sum(array2)
    
    return float(idx_2_only)/float(idx_2_total)


def critical_success(array1, array2):
    """
    calculate the critical success rate based upon 2 boolean maps. 
    """
    idx_both = np.sum(np.logical_and(array1, array2))
    idx_either = np.sum(np.logical_or(array1, array2))
    return float(idx_both)/float(idx_either)

def plot_contingency(x, y, contingency, title):
    """
    Prepare a geographical map of a contingency score map, with appropriate coloring
    """
    import cartopy.crs as ccrs
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import cartopy
    from matplotlib import colors
    from cartopy.io.img_tiles import Stamen
    import matplotlib.pyplot as plt
    plot_image = np.ma.masked_where(contingency==0, contingency)
    extent = (x.min(), x.max(), y.min(), y.max())
    cmap = colors.ListedColormap(['blue', 'red', 'green'])
    bounds=[0.5, 1.5, 2.5, 3.5]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    # get hold of the coastlines for that area.
#     ax.add_feature(cartopy.feature.LAND, zorder=1)
#     ax.add_feature(cartopy.feature.OCEAN, zorder=1)
#     ax.add_feature(cartopy.feature.COASTLINE)
#     ax.add_feature(cartopy.feature.BORDERS, linestyle=':')
#     ax.add_feature(cartopy.feature.LAKES, alpha=0.5)
#     ax.add_feature(cartopy.feature.RIVERS)
#    tiler = StamenTerrain()
    tiler = Stamen('terrain-background')
#    mercator = tiler.crs
    fig =plt.figure(figsize=(10,10))
#    ax = fig.add_subplot(1, 1, 1, projection=mercator)
    ax = fig.add_subplot(111, projection=ccrs.OSGB()) # mercator
#     ax.stock_img()

    

    ax.set_extent([231340, 245110, 829890, 842120], crs=ccrs.OSGB())
 #   ax.set_extent(extent)
    gl = ax.gridlines(crs=ccrs.OSGB(), draw_labels=False,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
#     ax.background_patch.set_fill(False)
#    ax.add_image(tiler, 7)
    img = ax.imshow(plot_image, extent=extent, vmin=1., vmax=3., interpolation='nearest', 
                        cmap=cmap, norm=norm, origin='lower', zorder=3, transform=ccrs.OSGB())  # origin='lower', transform=mercator
 #   ax.set_xlabel('longitude')
 #   ax.set_ylabel('latitude')


    # make a color bar
    cbar = plt.colorbar(img, cmap=cmap, norm=norm, boundaries=bounds, ticks=[1, 2, 3], orientation='horizontal')

    cbar.ax.set_xticklabels(['M1B0', 'M0B1', 'M1B1'])
    fig.suptitle(title, fontsize=14)
    fig.savefig('{:s}.png'.format(title), bbox_inches='tight', dpi=300)



def contingency(bench_fn, model_fn, bench_thres, model_thres, mask_fn, title, masking=False):
    #gis.gdal_warp(mask_fn, bench_fn, 'mask.tif')
    #gis.gdal_warp(model_fn, bench_fn, 'model.tif')
    #gis.gdal_warp(bench_fn, bench_fn, 'bench.tif')
    x, y, bench, fill_bench = gdal_readmap(bench_fn, 'GTiff')
    #x, y, model, fill_model = gdal_readmap('model.tif', 'GTiff')
    #x, y, mask, fill_mask = gdal_readmap('mask.tif', 'GTiff')
    x, y, model, fill_model = gdal_readmap(model_fn, 'GTiff')
    x, y, mask, fill_mask = gdal_readmap(mask_fn, 'GTiff')
    if masking:
        bench = np.ma.masked_where(np.logical_or(bench==fill_bench, mask==0), bench)
        model = np.ma.masked_where(np.logical_or(model==fill_model, mask==0), model)
    else:
        bench = np.ma.masked_where(bench==fill_bench, bench)
        model = np.ma.masked_where(model==fill_model, model)

    bench[bench==fill_bench] = 0.

    bench[model==fill_model] = 0.

    model[model==fill_model] = 0.
        
    flood1, flood2, cont_arr = contingency_map(bench, model, threshold1=bench_thres, threshold2=model_thres)
    hr = hit_rate(flood1, flood2)
    far = false_alarm_rate(flood1, flood2)
    csi = critical_success(flood1, flood2)

    return hr, far, csi, x, y, cont_arr

model_fn = "carlisle-5m.wd"

bench_fn = "carlisle-5m.dat"

mask_fn  = "mask.wd"

title = "Carlisle"

hr, far, csi, x, y, cont_arr = contingency(bench_fn, model_fn, 0.1, 0.1, mask_fn, title, masking = False)

print(y)
print(x)

print('Scores without urban mask')
print('Hit rate: {:f}'.format(hr))
print('False Alarm rate: {:f}'.format(far))
print('Critical success index: {:f}'.format(csi))


#xmin, xmax, ymin, ymax = [231340, 245110, 829890, 842120]
plot_contingency(x, y, np.flipud(cont_arr), 'Carlisle')
