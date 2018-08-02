import pyflann
from  netCDF4 import Dataset
import matplotlib.pyplot as plt
from skimage import measure
import math
import numpy as np
from scipy import spatial,io
import os
import timeit
from mpl_toolkits.basemap import Basemap
import inject_bathymetry
import mesh_definition_tools as mdt
import pprint
plt.switch_backend('agg')

######################################################################################
# Initial declarations
######################################################################################

km = 1000
deg2rad = np.pi/180
rad2deg = 180/np.pi

nn_search = "flann"
plot_option = True

# Region boxes
Delaware_Bay =  {"include":[np.array([-75.61903,-74.22, 37.767484, 40.312747])],
                 "exclude":[]}
Galveston_Bay = {"include":[np.array([-95.45,-94.4, 29, 30])],
                 "exclude":[]}
US_East_Coast = {"include":[np.array([-81.7,-62.3,25.1,46.24])],
                 "exclude":[np.array([-66.0,-64.0,31.5,33.0]),     # Bermuda
                             np.array([-79.75,-70.0,20.0,28.5])]}  # Bahamas
US_Gulf_Coast = {"include":[np.array([-98.0,-80.0,26.0,31.0])],
                 "exclude":[]}

US_West_Coast = {"include":[np.array([-127.0,-116.0,32.5,49.0])],
                 "exclude":[np.array([-116.5,-115.0,32.8,33.8]),  # Salton Sea
                            np.array([-120.5,-116.5,35.5,40.5])]} # Lake Tahoe, etc.

CONUS = {"include":[],"exclude":[]}
CONUS["include"].extend(US_East_Coast["include"]) 
CONUS["include"].extend(US_Gulf_Coast["include"]) 
CONUS["include"].extend(US_West_Coast["include"])
CONUS["exclude"].extend(US_East_Coast["exclude"])
CONUS["exclude"].extend(US_West_Coast["exclude"])
pprint.pprint(CONUS)


# Plotting boxes
Western_Atlantic = np.array([-98.186645, -59.832744, 7.791301 ,45.942453])
Contiguous_US = np.array([-132.0,-59.832744,7.791301,51.0])
Delware_Region = np.array([-77, -69.8 ,35.5, 41])
Entire_Globe = np.array([-180,180,-90,90])

######################################################################################
# User-defined inputs
######################################################################################

# Path to bathymetry data and name of file
data_path = "/users/sbrus/climate/bathy_data/SRTM15_plus/"
nc_file = "earth_relief_15s.nc"

# Bounding box of coastal refinement region
region_box = CONUS
origin = np.array([-100,40])

# Mesh parameters
grd_box = Entire_Globe 
ddeg = .1
dx_min = 10*km
dx_max = 60*km
trans_width = 600*km
trans_start = 400*km
mesh_type = 'EC'     #'EC' (defaults to 60to30), 'QU' (uses dx_max)

# Bounding box of plotting region
plot_box = Contiguous_US 

######################################################################################
# Functions
######################################################################################

def CPP_projection(lon,lat,origin):

  R = 6378206.4
  origin = origin*deg2rad
  x = R*(lon*deg2rad-origin[0])*np.cos(origin[1])
  y = R*lat*deg2rad

  return x,y

##############################################################

def create_background_mesh(grd_box,ddeg,mesh_type):

  # Create cell width background grid
  lat_grd = np.arange(grd_box[2],grd_box[3],ddeg)
  lon_grd = np.arange(grd_box[0],grd_box[1],ddeg)
  nx_grd = lon_grd.size
  ny_grd = lat_grd.size
  print "Background grid dimensions:", ny_grd,nx_grd
  
  # Assign background grid cell width values
  if mesh_type == 'QU':
    cell_width = dx_max*np.ones(D.shape)
  elif mesh_type == 'EC':
    cell_width_lat = mdt.EC_CellWidthVsLat(lat_grd)
    cell_width = np.tile(cell_width_lat,(nx_grd,1)).T*km
  
  # Plot background cell width
  if plot_option:
    plt.figure()
    plt.contourf(lon_grd,lat_grd,cell_width)
    m.drawcoastlines()
    plt.colorbar()
    plt.savefig('bckgnd_grid_cell_width.png',bbox_inches='tight')

  return (lon_grd,lat_grd,cell_width)

##############################################################

def get_data_inside_box(lon,lat,data,box):

  # Find indicies of coordinates inside bounding box
  lon_idx, = np.where((lon > box[0]) & (lon < box[1]))
  lat_idx, = np.where((lat > box[2]) & (lat < box[3]))
  
  # Get region data inside bounding box
  lon_region = lon[lon_idx]
  lat_region = lat[lat_idx]
  latlon_idx = np.ix_(lat_idx,lon_idx)
  dtype =  type(data)
  try:     # Numpy indexing
    z_region = data[latlon_idx]
  except:  # NetCDF indexing
    z_region = data[lat_idx,lon_idx]

  return (lon_region,lat_region,z_region)

##############################################################

def extract_coastlines(nc_file,region_box,z_contour=0,n_longest=10):

  # Open NetCDF file and read cooordintes
  nc_fid = Dataset(nc_file,"r")
  lon = nc_fid.variables['lon'][:]
  lat = nc_fid.variables['lat'][:]
  zdata = nc_fid.variables['z']
  
  # Get coastlines for refined region
  coastline_list = []
  for box in region_box["include"]:

    # Find coordinates and data inside bounding box
    lon_region,lat_region,z_region = get_data_inside_box(lon,lat,zdata,box)
    print "Regional bathymetry data shape:", z_region.shape
    
    # Find coastline contours
    print "Extracting coastline"
    contours = measure.find_contours(z_region,z_contour)

    # Keep only n_longest coastlines and those not within exclude areas
    contours.sort(key=len,reverse=True)
    for c in contours[:n_longest]:
      # Convert from pixel to lon,lat
      c[:,0] = (box[3]-box[2])/float(len(lat_region))*c[:,0] + box[2]        
      c[:,1] = (box[1]-box[0])/float(len(lon_region))*c[:,1] + box[0]
      c = np.fliplr(c)
  
      exclude = False
      for area in region_box["exclude"]:
        # Determine coastline coordinates in exclude area
        x_idx = np.where((c[:,0] > area[0]) & (c[:,0] < area[1]))
        y_idx = np.where((c[:,1] > area[2]) & (c[:,1] < area[3]))
        idx = np.intersect1d(x_idx,y_idx)
  
        # Exlude coastlines that are entirely contained in exclude area
        if idx.size == c.shape[0]:
          exclude = True
          break
  
      # Keep coastlines not entirely contained in exclude areas
      if not exclude:
        cpad= np.vstack((c,[np.nan,np.nan]))
        coastline_list.append(cpad)

    print "Done"
  
  # Combine coastlines 
  coastlines = np.concatenate(coastline_list)
  
  if plot_option:
  
    # Find coordinates and data inside plotting box
    lon_plot,lat_plot,z_plot = get_data_inside_box(lon,lat,zdata,plot_box)
 
    # Plot bathymetry data and coastlines
    plt.figure()
    levels = np.linspace(np.amin(z_plot),np.amax(z_plot),100)
    ds = 10                              # Downsample
    dsx = np.arange(0,lon_plot.size,ds)  # bathy data
    dsy = np.arange(0,lat_plot.size,ds)  # to speed up
    dsxy = np.ix_(dsy,dsx)               # plotting
    plt.contourf(lon_plot[dsx],lat_plot[dsy],z_plot[dsxy],levels=levels)
    m.drawcoastlines()
    plt.plot(coastlines[:,0],coastlines[:,1],color='white')
    plt.colorbar()
    plt.axis('equal')
    plt.savefig('bathy_coastlines.png',bbox_inches='tight')

  return coastlines

##############################################################

def distance_to_coast(coastlines,lon_grd,lat_grd,origin):

  # Convert to x,y and create kd-tree
  coast_pts = coastlines[np.isfinite(coastlines).all(axis=1)]
  coast_pts_xy = np.copy(coast_pts)
  coast_pts_xy[:,0],coast_pts_xy[:,1] = CPP_projection(coast_pts[:,0],coast_pts[:,1],origin)
  if nn_search == "kdtree":
    tree = spatial.KDTree(coast_pts_xy)
  elif nn_search == "flann":
    flann = pyflann.FLANN()
    flann.build_index(coast_pts_xy,algorithm='kdtree',target_precision=.999999999)
  
  # Put x,y backgound grid coordinates in a nx_grd x 2 array for kd-tree query
  Lon_grd,Lat_grd = np.meshgrid(lon_grd,lat_grd)
  X_grd,Y_grd = CPP_projection(Lon_grd,Lat_grd,origin)
  pts = np.vstack([X_grd.ravel(), Y_grd.ravel()]).T
  
  # Find distances of background grid coordinates to the coast
  print "Finding distance"
  start = timeit.default_timer()
  if nn_search == "kdtree":
    d,idx = tree.query(pts)
  elif nn_search == "flann":
    idx,d = flann.nn_index(pts,checks=500)
    d = np.sqrt(d)
  end = timeit.default_timer()
  print "Done"
  print end-start, " seconds"
  
  # Make distance array that corresponds with cell_width array
  D = np.reshape(d,Lon_grd.shape)
  
  if plot_option:

    # Find coordinates and data inside plotting box
    lon_plot,lat_plot,D_plot = get_data_inside_box(lon_grd,lat_grd,D,plot_box)
 
    # Plot distance to coast
    plt.figure()
    D_plot = D_plot/km 
    levels = np.linspace(np.amin(D_plot),np.amax(D_plot),10)
    plt.contourf(lon_plot,lat_plot,D_plot,levels=levels)
    m.drawcoastlines()
    plt.plot(coastlines[:,0],coastlines[:,1],color='white')
    plt.grid(xdata=lon_plot,ydata=lat_plot,c='k',ls='-',lw=0.1,alpha=0.5)
    plt.colorbar()
    plt.axis('equal')
    plt.savefig('distance.png',bbox_inches='tight')

  return D

##############################################################

def compute_cell_width(D,cell_width):

  # Compute cell width based on distance
  backgnd_weight = .5*(np.tanh((D-trans_start-.5*trans_width)/(.2*trans_width))+1)
  dist_weight = 1-backgnd_weight
  ## Use later for depth and slope dependent resolution
  ##hres = np.maximum(dx_min*bathy_grd/20,dx_min)
  ##hres = np.minimum(hres,dx_max)
  #hw = np.zeros(Lon_grd.shape) + dx_max
  #hw[ocn_idx] = np.sqrt(9.81*bathy_grd[ocn_idx])*12.42*3600/25
  #hs = .20*1/dbathy_grd
  #h = np.fmin(hs,hw)
  #h = np.fmin(h,dx_max)
  #h = np.fmax(dx_min,h)
  cell_width = (dx_min*dist_weight+np.multiply(cell_width,backgnd_weight))/km
  
  # Save matfile
  io.savemat('cellWidthVsLatLon.mat',mdict={'cellWidth':cell_width,'lon':lon_grd,'lat':lat_grd})
  
  if plot_option:

    # Find coordinates and data inside plotting box
    lon_plot,lat_plot,cell_width_plot = get_data_inside_box(lon_grd,lat_grd,cell_width,plot_box)

    # Plot cell width
    plt.figure()
    levels = np.linspace(np.amin(cell_width_plot),np.amax(cell_width_plot),100)
    plt.contourf(lon_plot,lat_plot,cell_width_plot,levels=levels)
    m.drawcoastlines()
    plt.plot(coastlines[:,0],coastlines[:,1],color='white')
    plt.colorbar()
    plt.axis('equal')
    plt.savefig('cell_width.png',bbox_inches='tight')

######################################################################################
# Main 
######################################################################################

if __name__ == "__main__":

  # Get coarse coastlines for ploting
  m = Basemap(projection='cyl',llcrnrlat=plot_box[2],urcrnrlat=plot_box[3],\
              llcrnrlon=plot_box[0],urcrnrlon=plot_box[1],resolution='c')
 
  # Create the background cell width array 
  lon_grd,lat_grd,cell_width = create_background_mesh(grd_box,ddeg,mesh_type)
 
  # Get coastlines from bathy/topo  data set   
  coastlines = extract_coastlines(data_path+nc_file,region_box)
  
  # Compute distance from background grid points to coastline 
  D = distance_to_coast(coastlines,lon_grd,lat_grd,origin)
  
  # Blend coastline and background resolution, save cell_width array as .mat file  
  compute_cell_width(D,cell_width)


## Incorporate later for depth and slope dependent resolution

## Interpolate bathymetry onto background grid
#Lon_grd = Lon_grd*deg2rad
#Lat_grd = Lat_grd*deg2rad
#bathy = inject_bathymetry.interpolate_bathymetry(data_path+nc_file,Lon_grd.ravel(),Lat_grd.ravel())
#bathy_grd = -1.0*np.reshape(bathy,(ny_grd,nx_grd))
#ocn_idx = np.where(bathy_grd > 0)
#
#if plot_option:
#  plt.figure()
#  levels = np.linspace(0,11000,100)
#  plt.contourf(lon_grd,lat_grd,bathy_grd,levels=levels)
#  m.drawcoastlines()
#  plt.colorbar()
#  plt.axis('equal')
#  plt.savefig('bckgnd_grid_bathy.png',bbox_inches='tight')

## Interpolate bathymetry gradient onto background grid
#dbathy = inject_bathymetry.interpolate_bathymetry(data_path+nc_file,Lon_grd.ravel(),Lat_grd.ravel(),grad=True)
#dbathy = np.reshape(dbathy,(ny_grd,nx_grd))
#dbathy_grd = np.zeros((ny_grd,nx_grd))
#dbathy_grd[ocn_idx] = dbathy[ocn_idx]
#
#if plot_option:
#  plt.figure()
#  plt.contourf(lon_grd,lat_grd,1/dbathy_grd)
#  m.drawcoastlines()
#  plt.colorbar()
#  plt.axis('equal')
#  plt.savefig('bckgnd_grid_bathy_grad.png',bbox_inches='tight')




