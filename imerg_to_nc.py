import numpy as np
import h5py
import xarray as xr
import matplotlib.pyplot as plt
import sys
import re
import pandas as pd

# Load the HDF5 file from command-line argument
input_file = sys.argv[1]
data = h5py.File(input_file, 'r')

# Extract the precipitation array
precip = data['/Grid/precipitation'][:]

# Remove the dummy first dimension, transpose, and flip vertically
precip = np.flip(precip[0, :, :].transpose(), axis=0)

# Define latitude and longitude arrays for plotting and storing in NetCDF
lons = np.linspace(-180, 180, precip.shape[1])
lats = np.linspace(-90, 90, precip.shape[0])

# Extract date and start time from filename using regex
match = re.search(r"\.(\d{8})-S(\d{6})", input_file)
if match:
    date_str = match.group(1)  # '20240409'
    time_str = match.group(2)  # '210000'

    # Convert to pandas datetime
    timestamp = pd.to_datetime(date_str + time_str, format='%Y%m%d%H%M%S')
    print("Converted time:", timestamp)
else:
    print("Date and time information not found in the filename.")
    sys.exit(1)

# Expand precipitation data to add a new time dimension
precip = precip[np.newaxis, :, :]  # Add a new axis for time, making it (1, latitude, longitude)


# Create an xarray DataArray and Dataset for the precipitation data with time
precip_da = xr.DataArray(
    precip,
    dims=["time", "latitude", "longitude"],
    coords={"time": [timestamp], "latitude": lats, "longitude": lons},
    attrs={"units": "mm/hr", "long_name": "precipitation rate"},
)

# Wrap into a Dataset and add global attributes
ds = xr.Dataset({"precipitation": precip_da})
ds.attrs["title"] = "IMERG Precipitation Data"
ds.attrs["source"] = "NASA IMERG"
ds.attrs["description"] = "Precipitation rate in mm/hr"
ds.attrs["time_coverage_start"] = str(timestamp)

# Define output NetCDF filename
output_file = "imerg_precipitation_with_time.nc"
output_file = input_file + ".nc"

# Write the Dataset to a NetCDF file
ds.to_netcdf(output_file)
print(f"Data written to {output_file}")
