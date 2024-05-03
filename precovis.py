#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd
import sqlite3 as sql
from astropy.time import Time
from astroquery.jplsbdb import SBDB

from precovery.orbit import Orbit
from precovery.orbit import EpochTimescale

import numpy as np
import healpy as hp

from datetime import date

from bokeh.plotting import figure, output_file, show, save
from bokeh.models import Range1d, Column, Row, CustomJS, DateRangeSlider, DateSlider, RangeSlider, Div, Slider, ColumnDataSource
from bokeh.models.tools import BoxZoomTool, ResetTool, HoverTool, TapTool
from bokeh.events import DoubleTap, ButtonClick
from bokeh.layouts import layout
from bokeh.io import curdoc
from bokeh.document import Document
from bokeh.models.widgets import DataTable, DateFormatter, TableColumn, Button, AbstractButton

from bokeh.models import ColorBar, ColumnDataSource
from bokeh.palettes import Spectral6
from bokeh.plotting import figure, output_file, show
from bokeh.transform import linear_cmap


NSIDE = 32
NPIX = hp.nside2npix(NSIDE)
con = sql.connect("index.db")
frames = pd.read_sql("""SELECT * FROM frames""", con)


def prepareOrbit(asteroid):
    """Create and return an orbit object of a given asteroid from the SBDB database"""
    result = SBDB.query(asteroid, full_precision=True, phys=True)

    # Extract the epoch into a Time object
    epoch = Time(result["orbit"]["epoch"], scale="tdb", format="jd")
    epoch_tt_mjd = epoch.tt.mjd

    # Extract physical characteristics if they exit
    if "phys_par" in result:
        if "H" in result["phys_par"]:
            H = result["phys_par"]["H"]
        else:
            H = 20.0

        if "G" in result["phys_par"]:
            G = result["phys_par"]["G"]
        else:
            G = 0.15

    # Define a precovery orbit
    orbit = Orbit.keplerian(
        0,
        result["orbit"]["elements"]["a"].value,
        result["orbit"]["elements"]["e"],
        result["orbit"]["elements"]["i"].value,
        result["orbit"]["elements"]["om"].value,
        result["orbit"]["elements"]["w"].value,
        result["orbit"]["elements"]["ma"].value,
        epoch_tt_mjd,
        EpochTimescale.TT,
        H, 
        G,
    )
    
    return orbit



def asteroidPath(orbit,start,end):
    """Return dataframe containing the ra, dec and time of a given object's orbit within the given timeframe
    
    Keyword arguments:
    orbit -- string name of the object
    start -- integer mjd value of the start of the time range for the orbit's path
    end -- integer mjd value of the end of the time range for the orbit's path
    """
    year_list = [x for x in np.linspace(start,end,20*(end-start))]
    ephemeris_list = orbit.compute_ephemeris("I11", year_list)
    ra = []
    dec = []
    time = []
    for ephemeris in ephemeris_list:
        ra.append(ephemeris.ra)
        dec.append(ephemeris.dec)
        time.append(ephemeris.mjd)

    coords_df = pd.DataFrame({"ra": ra, "dec": dec, "time": time})
    coords_df['healpixel'] = coordToPixel(coords_df['dec'].values,coords_df['ra'].values)
    
    return coords_df



def pixelToCoord(vectors):
    """Return the ra and dec tuple from a given healpixel pixel"""
    theta, phi = hp.pixelfunc.vec2ang(vectors)
    ra = np.degrees(np.pi*2.-phi)
    dec = -np.degrees(theta-np.pi/2.)
    return ra, dec



def coordToPixel(dec,ra):
    """Return the healpixel tuple from a given ra and dec"""
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-dec+90.),np.radians(360.-ra))



def pixelToBoundaries(pixels):
    """Return a dataframe with ra and dec coordinates for each healpixel boundary given within a list"""
    pixel_boundaries = {}
    for pixel in pixels:
        hp_boundaries = hp.boundaries(nside=NSIDE,pix=pixel,step=1).T
        ra,dec = pixelToCoord(hp_boundaries)
        boundary_path = pd.DataFrame({"ra": ra, "dec": dec})
        boundary_path = pd.concat([boundary_path,boundary_path.iloc[:1]], ignore_index=True)
        pixel_boundaries[pixel] = boundary_path
    return pixel_boundaries



def findGraphTransitions(coords):
    """Return a list of indexes where the ra and dec of a line with cross from one edge of the graph to the other"""
    large_diff = [0]
    for i in range(len(coords)-1):
        if np.abs(coords.iloc[i+1]['ra']-coords.iloc[i]['ra']) > 350:
            large_diff.append(i)
    large_diff.append(len(coords)-1)
    return large_diff



def mjdSliderCreation(startDate,endDate):
    """Create and return bokeh range slider with mjd dates. 

    Keyword arguments:
    startDate -- integer mjd value of the start of the time range of the slider
    endDate -- integer mjd value of the end of the time range of the slider
    """
    #startDate = 58138 #1/20/2018 to 4/29/2018
    duration = endDate-startDate
    datesX = [x for x in range(startDate,endDate+1)]
    
    date_range_slider = RangeSlider(
        title="MJD Date ("+str(duration)+" days)", start=datesX[0], end=datesX[duration],
        value=(datesX[0], datesX[duration]), step=1, width=300)
    
    return date_range_slider




def tableCreation(dataframe, fields):
    """Create and return bokeh table from dataframe. 

    Keyword arguments:
    dataframe -- dataframe where the table's information is populated from
    fields -- list of columns from dataframe that will be used within the table
    """
    global table_source
    table_source = ColumnDataSource(data=dataframe)
    global table_ref_source
    table_ref_source = ColumnDataSource(data=dataframe)

    columns = [TableColumn(field=x, title=x) for x in fields]

    data_table = DataTable(source=table_source, columns=columns, editable=True, height=500, width=900,index_width=60)
    return data_table




def plotIntercectingHealpixels(plot, objectName, startDate, endDate, dur):    
    """Plot healpixels that intercept with an object's orbit and date within the time range (startDate, endDate) as bokeh patches on a given plot. 

    Keyword arguments:
    plot -- the bokeh figure where the orbit's path is plotted
    objectName -- string name of the object
    startDate -- integer mjd value of the start of the time range of the time column
    endDate -- integer mjd value of the end of the time range of the time column (inclusive)
    dur -- float value of how similar the healpixels and the frames mjd values can be
    """
    
    #Orbit Set-up
    coords_df, transitions_orbit = objectToCoordsDF(objectName, startDate, endDate, transitions=True)


    #Healpixel Set-up
    ztf_frames = frames[(frames["dataset_id"] == "ztf")]
    ztf_frames_sorted = ztf_frames.sort_values("exposure_mjd_mid")
    merged_healpixels = pd.merge_asof(ztf_frames_sorted, coords_df, left_on="exposure_mjd_mid", right_on="time", by="healpixel", tolerance=dur).dropna() #by="healpixel"

    healpixel_time_groups = merged_healpixels.groupby('healpixel')
    pixels = list(healpixel_time_groups.groups.keys())
    #time_lists = list(healpixel_time_groups['exposure_mjd_mid'].apply(list))

    healpixels_df = pd.DataFrame(data={'healpixel':pixels, 'times':list(healpixel_time_groups['exposure_mjd_mid'].apply(list)), 'id':list(healpixel_time_groups['id'].apply(list)), 'dataset_id':list(healpixel_time_groups['dataset_id'].apply(list)), 'exposure_id':list(healpixel_time_groups['exposure_id'].apply(list)), 'obscode':list(healpixel_time_groups['obscode'].apply(list))})
    pixels_list = healpixels_df['healpixel'].tolist()
    
    #Healpixel Plotting
    patches = plotHealpixels(plot, pixels_list)
    
    return patches, merged_healpixels




'''
def plotHealpixels(plot, pixels):
    """Plot given healpixels as bokeh patches on a given plot. 

    Keyword arguments:
    plot -- the bokeh figure where the orbit's path is plotted
    pixels -- list of healpixels to plot 
    """
    
    pixel_boundaries = pixelToBoundaries(pixels)
    
    #Healpixel Plotting
    patches = [] #Due to edge cases, there will be more patches than actual healpixels (two patches for one edge healpixel)
    for pixel in pixel_boundaries.keys():
            if ((pixel_boundaries.get(pixel)['ra'].max() - pixel_boundaries.get(pixel)['ra'].min()) > 300):
                left_patch = pixel_boundaries.get(pixel)
                right_patch = pixel_boundaries.get(pixel).copy(deep=True)
                for p in range(len(left_patch)):
                    if (left_patch['ra'][p] > 350):
                        left_patch['ra'][p] = left_patch['ra'][p] - 360
                    if (right_patch['ra'][p] < 10):
                        right_patch['ra'][p] = right_patch['ra'][p] + 360        
                patches.append(plot.patch(left_patch['ra'],left_patch['dec'], line_width=0.5, line_color='white', color='skyblue'))#,line_color='white', color='skyblue'
                patches.append(plot.patch(right_patch['ra'],right_patch['dec'], line_width=0.5, line_color='white', color='skyblue'))#,line_color='white', color='skyblue'
            else:
                patches.append(plot.patch(pixel_boundaries.get(pixel)['ra'],pixel_boundaries.get(pixel)['dec'], line_width=0.5, line_color='white', color='skyblue'))#,line_color='white', color='skyblue'

    return patches
'''

def plotHealpixels(plot, pixels):
    """Plot given healpixels as bokeh patches on a given plot. 

    Keyword arguments:
    plot -- the bokeh figure where the orbit's path is plotted
    pixels -- list of healpixels to plot 
    """
    
    pixel_boundaries = pixelToBoundaries(pixels)
    
    #Healpixel Plotting
    patches = [] #Due to edge cases, there will be more patches than actual healpixels (two patches for one edge healpixel)
    global pixel_sources
    pixel_sources = []

    for pixel in pixel_boundaries.keys():
            
        pixel_boundary_df = pixel_boundaries.get(pixel)
        
        if ((pixel_boundary_df['ra'].max() - pixel_boundary_df['ra'].min()) > 300):
            left_patch = pixel_boundary_df
            right_patch = pixel_boundary_df.copy(deep=True)
            for p in range(len(left_patch)):
                if (left_patch['ra'][p] > 350):
                    left_patch['ra'][p] = left_patch['ra'][p] - 360
                if (right_patch['ra'][p] < 10):
                    right_patch['ra'][p] = right_patch['ra'][p] + 360 
            left_pixel_source = ColumnDataSource(left_patch)
            right_pixel_source = ColumnDataSource(right_patch)
            pixel_sources.append(left_pixel_source)
            pixel_sources.append(right_pixel_source)
            patches.append(plot.patch('ra','dec', source=left_pixel_source, line_width=0.5, line_color='white', color='skyblue'))#,line_color='white', color='skyblue'
            patches.append(plot.patch('ra','dec', source=right_pixel_source, line_width=0.5, line_color='white', color='skyblue'))#,line_color='white', color='skyblue'
        else:
            pixel_source = ColumnDataSource(pixel_boundary_df)
            pixel_sources.append(pixel_source)
            patches.append(plot.patch('ra','dec', source=pixel_source, line_width=0.5, line_color='white', color='skyblue'))#,line_color='white', color='skyblue'

    tap_tool = TapTool(renderers=patches)
    plot.add_tools(tap_tool)
            
    return patches



def plotOrbitPaths(plot, objectName, startDate, endDate):
    """Plot the path of orbit of a given object on a given plot with RA and Dec.

    Keyword arguments:
    plot -- the bokeh figure where the orbit's path is plotted
    objectName -- string name of the object
    startDate -- integer mjd value of the start of the time range of the time column
    endDate -- integer mjd value of the end of the time range of the time column (inclusive)
    """
    
    #Orbit Set-up
    coords_df, transitions_orbit = objectToCoordsDF(objectName, startDate, endDate, transitions=True)

    #Orbit Plotting
    orbit_lines = []
    global sources
    sources = []
    global ref_sources
    ref_sources = []
    global all_sources
    all_sources = []
    #mapper = linear_cmap(field_name='time', palette=Spectral6, low=start, high=end)
    for i in range(len(transitions_orbit)-1):
        sources.append(ColumnDataSource(coords_df.iloc[transitions_orbit[i]+1:transitions_orbit[i+1]]))
        ref_sources.append(ColumnDataSource(coords_df.iloc[transitions_orbit[i]+1:transitions_orbit[i+1]].copy(deep=True)))
        all_sources.append(sources[i])
        all_sources.append(ref_sources[i])
        orbit_lines.append(plot.line(x='ra',y='dec', source=sources[i], legend_label="Orbit Path", color="steelblue", line_width=1.0))
        plot.add_tools(HoverTool(tooltips=[("Ra", "@ra"),("Dec", "@dec"),("Time", "@time{0.00}"),("Healpixel", "@healpixel")],mode = "mouse",renderers=[orbit_lines[i]]))
    return orbit_lines, all_sources




def objectToCoordsDF(objectName, startDate, endDate, transitions=False):
    """Return a dataframe with RA, Dec, healpixel, and time (mjd) columnns of an object's path from the times startDate to endDate.

    Keyword arguments:
    objectName -- string name of the object
    startDate -- integer mjd value of the start of the time range of the time column
    endDate -- integer mjd value of the end of the time range of the time column (inclusive)
    transitions -- if True, return both the dataframe and a list of indeces where the RA goes between 0 and 360 between two rows of the dataframe
    """
    orbit = prepareOrbit(objectName) #Possible orbits: "IVEZIC" or "2020 AV2"
    coords_df = asteroidPath(orbit,startDate,endDate)
    
    if transitions:
        return coords_df, findGraphTransitions(coords_df)
    return coords_df




def plotSetUp(title):
    """Create a plot with with (0,360) RA as x-axis and (-90,90) Dec as y-axis named the variable title."""
    plot = figure(title=title, x_axis_label='ra', y_axis_label='dec', tools=[BoxZoomTool(), ResetTool()])
    plot.x_range = Range1d(0, 360)
    plot.y_range = Range1d(-90, 90)
    return plot

