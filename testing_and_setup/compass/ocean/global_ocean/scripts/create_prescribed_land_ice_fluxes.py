#!/usr/bin/env python
"""
This script creates forcing files for prescribed ice-shelf melt fluxes inferred
from satellite observations
"""

import configparser
import numpy as np
import xarray as xr
import sys
from mpas_tools.io import write_netcdf
import pyproj
from pyremap import MpasMeshDescriptor, ProjectionGridDescriptor, Remapper


def main():

    config = configparser.ConfigParser(
        interpolation=configparser.ExtendedInterpolation())
    config.read("config_prescribed_land_ice_fluxes.ini")
    prescribed_land_ice_fluxes(config)


def prescribed_land_ice_fluxes(config):

    in_filename = config.get('main', 'rignot_2013_file')

    out_filename = 'prescriped_land_ice_fluxes_rignot2013.nc'

    mpiTasks = config.getint('main', 'nprocs')

    remap_rignot(inFileName=in_filename, meshFileName='mesh.nc',
                 meshName='mpas_mesh', outFileName=out_filename,
                 mappingDirectory='.', method='conserve',
                 renormalizationThreshold=None, inVarName='melt_actual',
                 mpiTasks=mpiTasks)


def remap_rignot(inFileName, meshFileName, meshName, outFileName,
                 mappingDirectory='.', method='conserve',
                 renormalizationThreshold=None, inVarName='melt_actual',
                 mpiTasks=1):
    """
    Remap the Rignot et al. (2013) melt rates at 1 km resolution to an MPAS
    mesh

    Parameters
    ----------
    inFileName : str
        The original Rignot et al. (2013) melt rates

    meshFileName : str
        The MPAS mesh

    meshName : str
        The name of the mesh (e.g. oEC60to30wISC), used in the name of the
        mapping file

    outFileName : str
        The melt rates interpolated to the MPAS mesh with ocean sensible heat
        fluxes added on (assuming insulating ice)

    mappingDirectory : str
        The directory where the mapping file should be stored (if it is to be
        computed) or where it already exists (if not)

    method : {'bilinear', 'neareststod', 'conserve'}, optional
        The method of interpolation used, see documentation for
        `ESMF_RegridWeightGen` for details.

    renormalizationThreshold : float, optional
        The minimum weight of a denstination cell after remapping, below
        which it is masked out, or ``None`` for no renormalization and
        masking.

    inVarName : {'melt_actual', 'melt_steadystate'}
        Whether to use the melt rate for the time period covered in Rignot et
        al. (2013) with observed thinning/thickening or the melt rates that
        would be required if ice shelves were in steady state.

    mpiTasks : int, optional
        The number of MPI tasks to use to compute the mapping file
    """

    ds = xr.open_dataset(inFileName)
    lx = np.abs(1e-3 * (ds.xaxis.values[-1] - ds.xaxis.values[0]))
    ly = np.abs(1e-3 * (ds.yaxis.values[-1] - ds.yaxis.values[0]))

    inGridName = '{}x{}km_1.0km_Antarctic_stereo'.format(lx, ly)

    projection = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 '
                             '+k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')

    inDescriptor = ProjectionGridDescriptor.read(
        projection,  inFileName, xVarName='xaxis', yVarName='yaxis',
        meshName=inGridName)

    # convert to the units and variable names expected in MPAS-O
    rho_fw = 1000.
    s_per_yr = 365.*24.*60.*60.
    latent_heat_of_fusion = 3.337e5
    freshwaterFlux = ds[inVarName]*rho_fw/s_per_yr
    freshwaterFlux = freshwaterFlux.expand_dims('Time', axis=0)
    ds['prescribedLandIceFreshwaterFlux'] = freshwaterFlux
    ds['prescribedLandIceHeatFlux'] = latent_heat_of_fusion * freshwaterFlux
    ds = ds.drop_vars(['melt_actual', 'melt_steadystate', 'lon', 'lat'])

    outDescriptor = MpasMeshDescriptor(meshFileName, meshName)

    mappingFileName = '{}/map_{}_to_{}.nc'.format(mappingDirectory, inGridName,
                                                  meshName)

    remapper = Remapper(inDescriptor, outDescriptor, mappingFileName)

    remapper.build_mapping_file(method=method, mpiTasks=mpiTasks,
                                tempdir='.')

    dsRemap = remapper.remap(
        ds, renormalizationThreshold=renormalizationThreshold)

    for field in ['prescribedLandIceFreshwaterFlux',
                  'prescribedLandIceHeatFlux']:
        # zero out the field where it's currently NaN
        dsRemap[field] = dsRemap[field].where(dsRemap[field].notnull(), 0.)

    dsRemap.attrs['history'] = ' '.join(sys.argv)
    write_netcdf(dsRemap, outFileName)


if __name__ == '__main__':
    main()
