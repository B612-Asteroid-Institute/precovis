## Precovis

Tools to visualize searches performed by [Iterative Precovery and Orbit Determination (IPOD)](https://github.com/B612-Asteroid-Institute/ipod) and [precovery](https://github.com/B612-Asteroid-Institute/precovery). 

Goals:  
1. Visualize an orbit's trajectory on the sky. Propagate the orbit in a series of time steps between user defined start and end times, and display the orbit's position as a function of time.  
2. Visualize the frames table as a function of time. Display healpix frames being observed as a function of survey time. To start, lets create a single visualization for each observatory code. 
3. Combine 1 and 2, and display the orbit's position as a function of time, and overlay the frames table on top of the orbit highlighting the frames where a potential intersection may have occured.
4. Create a multi-view dashboard that allows the user to interact the above visualizations.
5. Create a visualization that displays the individual point-source measurements within each frame and shows summary statistics such as the astrometric position and uncertainty (if available), measurement time, magnitude and filter. Add this visualization to the dashboard.
6. Devise a visualization to incorporate information about the residuals between the orbit's predicted position and the observations within each frame. Add this visualization to the dashboard.
7. Stretch : Reference cutouts for each observation on disk and visualize them in the dashboard when clicking on a point-source measurement.
8. Stretch : Select observations to be used in the orbit determination process and produce a list of the observations.
9. Stretch : ...

Tools:  
`precovery` - A Python package to perform precovery searches on a set of observations  
`ipod` - A Python package to perform iterative precovery and orbit determination on a set of observations. IPOD performs iterative searches for new observations by mapping the orbit's approximate on-sky uncertainty to the observations. As new observations are found, the orbit is refined and the search is performed again.   
`adam_core` - A Python package that defines a common set of utilies underpinning Asteroid Institute's open source tools.  
`bokeh` - Current visualization package of choice for individual visualizations.

Unknowns:
Possible visualization tools for individual visualizations and multi-view dashboards.
- `bokeh` - Bokeh server allows for interactive visualizations to be created and served.
- `dash` - A Python framework for building analytical web applications. 
- `plotly` - A Python graphing library that makes interactive, publication-quality graphs online.
- `panel` - A high-level app and dashboarding solution for Python. 
- `d3.js` - A JavaScript library for producing dynamic, interactive data visualizations in web browsers.
- others?
