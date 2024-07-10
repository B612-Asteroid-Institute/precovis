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


### Examples

To run these examples, you'll need to install the appropriate packages. Ideally, in a fresh virtual environment (such as one created by `conda` or `virtualenv`).

```bash
conda create -n precovis_py311 python=3.11
conda activate precovis_py311
pip install -r requirements.txt
```

#### Introducing the Orbits Class (Using External Data)
```python
from adam_core.orbits.query import query_sbdb

object_ids = ["2013 RR165", "Eros", "Apophis"]
orbits = query_sbdb(object_ids)
```

This returns an `Orbits` object (which itself is a `quivr` table). To convert the table to a `pandas` dataframe, you can do the following:

```python 
orbits.to_dataframe()
```

`quivr` tables can be sliced like numpy arrays. 

```python
orbit0 = orbits[0]
orbit1 = orbits[1:2]
orbit2 = orbits[-1]
```

`quivr` tables also have several very useful functions. If you want to select the rows that match a particular value:
```python
filtered_orbits = orbits.select("object_id", "(2013 RR165)")
```

You can also apply masks. These are a bit more trickier since they require `pyarrow` compute functions. Here is an example:
```python
import pyarrow.compute as pc

filtered_orbits = orbits.apply_mask(pc.equal(orbits.object_id, "(2013 RR165)")) # This is equivalent to the select statement above
```

Here is a more complex one:
```python
filtered_orbits = orbits.apply_mask(
    pc.and_(
        pc.greater_equal(orbits.coordinates.x, 0.5), 
        pc.less_equal(orbits.coordinates.y, 0.0))
    )
```

Let's take a look at the actual definition of the `Orbits` object which we repeat here for convenience:

```python
import quivr as qv
from adam_core.coordinates import CartesianCoordinates

class Orbits(qv.Table):

    orbit_id = qv.LargeStringColumn(default=lambda: uuid.uuid4().hex)
    object_id = qv.LargeStringColumn(nullable=True)
    coordinates = CartesianCoordinates.as_column()
```

Notice how in adam_core, the default representation for all orbits is Cartesian. [CartesianCoordinates](#https://github.com/B612-Asteroid-Institute/adam_core/blob/main/src/adam_core/coordinates/cartesian.py#L20) are just another `quivr` table that is nested within the `Orbits` table. 

#### Defining An Orbit

Let's define an orbit without using external data. 

```python
from adam_core.orbits import Orbits
from adam_core.time import Timestamp
from adam_core.coordinates import KeplerianCoordinates, Origin

keplerian = KeplerianCoordinates.from_kwargs(
    a=[2.0],
    e=[0.1],
    i=[0.1],
    raan=[10.0],
    ap=[30.0],
    M=[23.0],
    time=Timestamp.from_kwargs(
        days=[59000],
        nanos=[0],
        scale="tdb"
    ),
    origin=Origin.from_kwargs(code=["SUN"]),
    frame="ecliptic",
)

orbits = Orbits.from_kwargs(
    object_id=["Dummy Object"],
    coordinates=keplerian.to_cartesian()
)
```

There is a lot going on here. Let's break it down:

1. We've defined a `KeplerianCoordinates` object. This object is just another `quivr` table. `KeplerianCoordinates` are a little easier to interpret geometrically than `CartesianCoordinates` so we use them this time to define an orbit. When we actually go to create the orbits class, we convert the `KeplerianCoordinates` to `CartesianCoordinates`.

2. The `.from_kwargs` constructor is one way to a create a `quivr` table from keywords representing the different columns of the underlying table. Notice how these keywords are lists, that is they represent columns and therefore are never scalars. Valid options are lists, `~numpy.ndarrays`, `~pandas.Series`, `~pyarrow.Array`, `~quivr.Table`, among others.

3. The `KeplerianCoordinates` object has an associated time. The default class that represents time is the [Timestamp]() class. This class is also a `quivr` table. We've chosen to use our own custom class because under the hood it stores time as two integers (days and nanoseconds). The Timestamp class has convenience functions to convert to and from `~astropy.time.Time` objects (`Timestamp.from_astropy` and `Timestamp.to_astropy`).

4. Coordinates also need to have a defined Origin and Frame. The Origin is a `quivr` table that has a single column `code` that represents the origin of the coordinates. The Frame is a string that represents the frame of the coordinates. The frame is an attribute and not a column because it is the same for all coordinates in the table (this is why we can pass a single value and not a list of values).

### Propagating Orbits

Now, we get to the fun parts of dealing with orbits -- propagation!

adam_core has a Propagator class that wraps different propagators. For now we will stick with PYOORB:
```python
from adam_core.propagator.adam_pyoorb import PYOORBPropagator
from adam_core.time import Timestamp
import numpy as np

propagator = PYOORBPropagator()
times = Timestamp.from_mjd(np.linspace(59000, 59100, 100))

propagated_orbits = propagator.propagate_orbits(orbits, times)
```

This returns an `Orbits` object that has M x N rows where M is the number of times and N is the number of orbits. Again, you can use `propagated_orbits.to_dataframe()` to convert to a pandas dataframe and inspect the results.

### Ephemeris Generation

We can also generate predicted ephemerides for the orbits. This is useful for visualizing the orbits on the sky. Here, we need to define an observer or set of observers and pass them to the propagator. `I41` is the Zwicky Transient Facility (ZTF) observatory code. To define multiple observers (different observatory codes with different times), you would use `quivr.concatenate` to concatenate the observers into a single `Observers`.

```python
from adam_core.propagator.adam_pyoorb import PYOORBPropagator
from adam_core.time import Timestamp
from adam_core.observers import Observers
import numpy as np

propagator = PYOORBPropagator()
observers = Observers.from_code(
    "I41",
    Timestamp.from_mjd(np.linspace(59000, 59100, 100)),
)

ephemeris = propagator.generate_ephemeris(propagated_orbits, observers)
```
