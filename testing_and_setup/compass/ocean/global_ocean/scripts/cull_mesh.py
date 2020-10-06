#!/usr/bin/env python
"""
This script performs the first step of initializing the global ocean.  This
includes:
1. combining Natural Earth land coverage north of 60S with Antarctic
   ice coverage or grounded ice coverage from Bedmap2.
2. combining transects defining critical passages.*
3. combining points used to seed a flood fill of the global ocean.
4. create masks from land coverage.
5. add land-locked cells to land coverage mask.
6. create masks from transects.*
7. cull cells based on land coverage but with transects present
8. create flood-fill mask based on seeds
9. cull cells based on flood-fill mask
10. create masks from transects on the final culled mesh*
* skipped if flag --with_critical_passages not present

The optional --with_cavities flag indicates that ice-shelf cavities are present
and the grounded-ice mask from Bedmap2 should be used. The optional
--with_critical_passages flag indicates that critical passages are
to be opened. Otherwise, steps 2, 5 and 9 are skipped
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

from argparse import ArgumentParser
import xarray
import logging
import sys
import multiprocessing

from geometric_features import GeometricFeatures, FeatureCollection, \
    read_feature_collection
from mpas_tools.mesh import conversion
from mpas_tools.mesh.mask import compute_mpas_region_masks, \
    compute_mpas_transect_masks, compute_mpas_flood_fill_mask
from mpas_tools.cime.constants import constants
from mpas_tools.io import write_netcdf
from mpas_tools.ocean.coastline_alteration import widen_transect_edge_masks, \
    add_critical_land_blockages, add_land_locked_cells_to_mask
from mpas_tools.viz.paraview_extractor import extract_vtk


def compute_land_mask(gf, dsBaseMesh, logger, pool, with_cavities,
                      use_mesh_conversion_mask):

    # start with the land coverage from Natural Earth
    fcLandCoverage = gf.read(componentName='natural_earth', objectType='region',
                             featureNames=['Land Coverage'])

    # remove the region south of 60S so we can replace it based on ice-sheet
    # topography
    fcSouthMask = gf.read(componentName='ocean', objectType='region',
                          featureNames=['Global Ocean 90S to 60S'])

    fcLandCoverage = fcLandCoverage.difference(fcSouthMask)

    # Add "land" coverage from either the full ice sheet or just the grounded
    # part
    if with_cavities:
        fcAntarcticLand = gf.read(componentName='bedmachine',
                                  objectType='region',
                                  featureNames=['AntarcticGroundedIceCoverage'])
    else:
        fcAntarcticLand = gf.read(componentName='bedmachine',
                                  objectType='region',
                                  featureNames=['AntarcticIceCoverage'])

    fcLandCoverage.merge(fcAntarcticLand)

    # combine the two regions into a single feature for faster processing
    fcLandCoverage = fcLandCoverage.combine('Land Coverage')

    # save the feature collection to a geojson file
    fcLandCoverage.to_geojson('land_coverage.geojson')

    if use_mesh_conversion_mask:
        dsLandMask = conversion.mask(dsBaseMesh, fcMask=fcLandCoverage)
    else:
        dsLandMask = compute_mpas_region_masks(
            dsMesh=dsBaseMesh, fcMask=fcLandCoverage, maskTypes=('cell',),
            logger=logger, pool=pool, showProgress=True)

    dsLandMask = add_land_locked_cells_to_mask(dsLandMask, dsBaseMesh,
                                               latitude_threshold=43.0,
                                               nSweeps=20)

    return dsLandMask


def add_land_blockages(gf, dsBaseMesh, dsLandMask, earth_radius, logger, pool,
                       custom_land_blockages, with_critical_passages,
                       use_mesh_conversion_mask):

    land_blockages = with_critical_passages or \
        (custom_land_blockages is not None)

    if land_blockages:
        if with_critical_passages:
            # merge transects for critical land blockages into
            # critical_land_blockages.geojson
            fcCritBlockages = gf.read(componentName='ocean',
                                      objectType='transect',
                                      tags=['Critical_Land_Blockage'])
        else:
            fcCritBlockages = FeatureCollection()

        if custom_land_blockages is not None:
            fcCritBlockages.merge(read_feature_collection(
                custom_land_blockages))

        # create masks from the transects
        if use_mesh_conversion_mask:
            dsCritBlockMask = conversion.mask(dsBaseMesh,
                                              fcMask=fcCritBlockages)
        else:
            fcCritBlockages = fcCritBlockages.combine('Critical Land Blockages')
            dsCritBlockMask = compute_mpas_transect_masks(
                dsMesh=dsBaseMesh, fcMask=fcCritBlockages,
                earthRadius=earth_radius, maskTypes=('cell',), logger=logger,
                pool=pool, showProgress=True)
        dsLandMask = add_critical_land_blockages(dsLandMask, dsCritBlockMask)

    return dsLandMask


def get_critical_passages(gf, custom_critical_passages, with_critical_passages,
                          use_mesh_conversion_mask):

    critical_passages = with_critical_passages or \
        (custom_critical_passages is not None)

    if critical_passages:
        fcCritPassages = FeatureCollection()
        if with_critical_passages:
            # merge transects for critical passages into
            # critical_passages.geojson
            fcCritPassages.merge(gf.read(componentName='ocean',
                                         objectType='transect',
                                         tags=['Critical_Passage']))

        if custom_critical_passages is not None:
            fcCritPassages.merge(read_feature_collection(
                custom_critical_passages))

        # create masks from the transects
        if not use_mesh_conversion_mask:
            fcCritPassages = fcCritPassages.combine('Critical Passages')
    else:
        fcCritPassages = None

    return fcCritPassages


def add_critical_passages(dsMesh, fcCritPassages, earth_radius, logger, pool,
                          use_mesh_conversion_mask):

    if fcCritPassages is None:
        dsCritPassMask = None
    else:
        # create masks from the transects
        if use_mesh_conversion_mask:
            dsCritPassMask = conversion.mask(dsMesh, fcMask=fcCritPassages)
        else:
            fcCritPassages = fcCritPassages.combine('Critical Passages')
            dsCritPassMask = compute_mpas_transect_masks(
                dsMesh=dsMesh, fcMask=fcCritPassages,
                earthRadius=earth_radius, maskTypes=('cell', 'edge'),
                logger=logger, pool=pool, showProgress=True)

        # Alter critical passages to be at least two cells wide, to avoid sea ice
        # blockage.
        dsCritPassMask = widen_transect_edge_masks(dsCritPassMask, dsMesh,
                                                   latitude_threshold=43.0)

    return dsCritPassMask


def do_flood_fill(gf, dsCulledMesh, logger, use_mesh_conversion_mask):

    # create seed points for a flood fill of the ocean
    # use all points in the ocean directory, on the assumption that they are, in
    # fact, in the ocean
    fcSeed = gf.read(componentName='ocean', objectType='point',
                     tags=['seed_point'])

    # create a mask for the flood fill seed points
    if use_mesh_conversion_mask:
        dsSeedMask = conversion.mask(dsCulledMesh, fcSeed=fcSeed)
    else:
        dsSeedMask = compute_mpas_flood_fill_mask(dsCulledMesh, fcSeed=fcSeed,
                                                  logger=logger)

    # cull the mesh a second time using a flood fill from the seed points
    dsCulledMesh = conversion.cull(dsCulledMesh, dsInverse=dsSeedMask,
                                   graphInfoFileName='culled_graph.info')
    return dsCulledMesh


def create_culled_mesh(gf, logger, pool, earth_radius, netcdfFormat,
                       with_cavities, custom_land_blockages,
                       custom_critical_passages, with_critical_passages,
                       use_mesh_conversion_mask, preserve_floodplain):
    # Create the land mask based on the land coverage, i.e. coastline data.
    dsBaseMesh = xarray.open_dataset('base_mesh.nc')

    dsLandMask = compute_land_mask(gf, dsBaseMesh, logger, pool,
                                   with_cavities,
                                   use_mesh_conversion_mask)

    dsLandMask = add_land_blockages(gf, dsBaseMesh, dsLandMask, earth_radius,
                                    logger, pool, custom_land_blockages,
                                    with_critical_passages,
                                    use_mesh_conversion_mask)

    fcCritPassages = get_critical_passages(gf, custom_critical_passages,
                                           with_critical_passages,
                                           use_mesh_conversion_mask)

    dsCritPassMask = add_critical_passages(dsBaseMesh, fcCritPassages,
                                           earth_radius, logger, pool,
                                           use_mesh_conversion_mask)

    dsPreserve = []
    if dsCritPassMask is not None:
        dsPreserve.append(dsCritPassMask)

    if preserve_floodplain:
        dsPreserve.append(dsBaseMesh)

    # cull the mesh based on the land mask
    dsCulledMesh = conversion.cull(dsBaseMesh, dsMask=dsLandMask,
                                   dsPreserve=dsPreserve)

    dsCulledMesh = do_flood_fill(gf, dsCulledMesh, logger,
                                 use_mesh_conversion_mask)

    write_netcdf(dsCulledMesh, 'culled_mesh.nc', format=netcdfFormat)

    return dsCulledMesh, fcCritPassages


def write_ice_mask(gf, dsCulledMesh, logger, pool, netcdfFormat, with_cavities,
                   use_mesh_conversion_mask):
    if with_cavities:
        fcAntarcticIce = gf.read(componentName='bedmachine', objectType='region',
                                 featureNames=['AntarcticIceCoverage'])
        fcAntarcticIce.to_geojson('ice_coverage.geojson')
        if use_mesh_conversion_mask:
            dsMask = conversion.mask(dsCulledMesh, fcMask=fcAntarcticIce)
        else:
            dsMask = compute_mpas_region_masks(
                dsCulledMesh, fcMask=fcAntarcticIce, maskTypes=('cell',),
                logger=logger, pool=pool, showProgress=True)
        landIceMask = dsMask.regionCellMasks.isel(nRegions=0)
        dsLandIceMask = xarray.Dataset()
        dsLandIceMask['landIceMask'] = landIceMask

        write_netcdf(dsLandIceMask, 'land_ice_mask.nc', format=netcdfFormat)

        dsLandIceCulledMesh = conversion.cull(dsCulledMesh, dsMask=dsMask)
        write_netcdf(dsLandIceCulledMesh, 'no_ISC_culled_mesh.nc',
                     format=netcdfFormat)


def main():
    parser = ArgumentParser()
    parser.add_argument("--with_cavities", action="store_true",
                        dest="with_cavities",
                        help="Whether the mesh should include Antarctic "
                             "ice-shelf cavities")
    parser.add_argument("--with_critical_passages", action="store_true",
                        dest="with_critical_passages",
                        help="Whether the mesh should open the standard "
                             "critical passages and close land blockages from "
                             "geometric_features")
    parser.add_argument(
        "--custom_critical_passages", dest="custom_critical_passages",
        required=False,
        help="A geojson file with critical passages to open.  This "
        "file may be supplied in addition to or instead of "
        "the default passages (--with_critical_passages)")
    parser.add_argument(
        "--custom_land_blockages", dest="custom_land_blockages",
        required=False,
        help="A geojson file with critical land blockages to close. "
        "This file may be supplied in addition to or instead of "
        "the default blockages (--with_critical_passages)")
    parser.add_argument("--preserve_floodplain", action="store_true",
                        dest="preserve_floodplain",
                        help="Whether to use the cellSeedMask field in the "
                             "base mesh to preserve a floodplain at elevations "
                             "above z=0")
    parser.add_argument("--use_mesh_conversion_mask", action="store_true",
                        dest="use_mesh_conversion_mask",
                        help="Whether to use mpas_tools.mesh.conversion.mask() "
                             "and therefore the c++ mask creator rather than "
                             "the python-based mask creator")
    parser.add_argument("--mask_creation_process_count", required=False,
                        dest="mask_creation_process_count", type=int,
                        help="The number of processes to use to compute masks. "
                             "The default is to use all available cores")
    args = parser.parse_args()

    pool = None
    if not args.use_mesh_conversion_mask:
        multiprocessing.set_start_method('fork')
        process_count = args.mask_creation_process_count
        if process_count is None:
            process_count = multiprocessing.cpu_count()
        else:
            process_count = min(process_count, multiprocessing.cpu_count())

        if process_count > 1:
            pool = multiprocessing.Pool(process_count)

    # required for compatibility with MPAS
    netcdfFormat = 'NETCDF3_64BIT'

    earth_radius = constants['SHR_CONST_REARTH']

    # set up a logger for masking functions that don't produce output without
    # one
    logger = logging.getLogger(__name__)
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    logger.propagate = False

    critical_passages = args.with_critical_passages or \
        (args.custom_critical_passages is not None)

    gf = GeometricFeatures()

    dsCulledMesh, fcCritPassages = create_culled_mesh(
        gf, logger, pool, earth_radius, netcdfFormat, args.with_cavities,
        args.custom_land_blockages, args.custom_critical_passages,
        args.with_critical_passages, args.use_mesh_conversion_mask,
        args.preserve_floodplain)

    if fcCritPassages is not None:
        dsCritPassMask = add_critical_passages(dsCulledMesh, fcCritPassages,
                                               earth_radius, logger, pool,
                                               args.use_mesh_conversion_mask)
        write_netcdf(dsCritPassMask, 'critical_passages_mask_final.nc',
                     format=netcdfFormat)

    write_ice_mask(gf, dsCulledMesh, logger, pool, netcdfFormat,
                   args.with_cavities, args.use_mesh_conversion_mask)

    extract_vtk(ignore_time=True, dimension_list=['maxEdges='],
                variable_list=['allOnCells'], filename_pattern='culled_mesh.nc',
                out_dir='culled_mesh_vtk')

    if args.with_cavities:
        extract_vtk(ignore_time=True, dimension_list=['maxEdges='],
                    variable_list=['allOnCells'],
                    filename_pattern='no_ISC_culled_mesh.nc',
                    out_dir='no_ISC_culled_mesh_vtk')

    if pool is not None:
        pool.terminate()


if __name__ == '__main__':
    main()
