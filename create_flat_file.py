"""Create flattened shapefile of non-overlapping polygons.

This script creates a flat file, in which there are no overlapping polygons. This
is useful for calculating stats where only visible parts of many overlapping polygons
must be considered. For example, a flat file is useful in calculating sums of visible
area by acquisition date in of a mosaic of images of many different acquisition dates.

This is achieved by rasterizing the input polygons (using a unique id field, in this
case, 'unique_id') while storing the attribute table in a dictionary, then
re-vectorizing the raster and re-joining the attributes. It is much faster than having
to break all the polygons, and with a high resolution the accuracy trade-off is
negligable.

Input shapefile features should be already sorted in the order of which they will be
flattened. 

Usage:
  python create_flat_file.py --input_shp /path/to/input_shp.shp --flat_output_shp
  /path/to/flat_output_shp.shp --unique_id_field id
"""

import datetime
import os
import sys
import tempfile
import logging
logging.getLogger().addHandler(logging.StreamHandler())

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

from argparse import ArgumentParser

_MAX_VALUE_32BIT = 2 ** 31 - 1
_SRS_EPSG = 4326  # Lat Lon WGS84.
_CYLINDRICAL_EQUAL_AREA_PROJ4_DEF = (
    '+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +datum='
    'WGS84 +units=m +no_defs')

class Error(Exception):

  """General error raises by this module."""


class FlatFile(object):

  """Flat file class."""

  def __init__(self, output_flat_file):
    """Initialization.

    Args:
      output_flat_file: Path to output flat shapefile. Flat file is stack
      ranked according to 'unique_id' (low to high) where underneath and
      overlapping areas are merged.
    """
    self._output_flat_file = output_flat_file
    self._driver = ogr.GetDriverByName('ESRI Shapefile')

  def FlatLayerFromRasterize(self, input_shp, output_tif, output_flat_shp,
                             width=70000, height=35000, attribute='unique_id',
                             x_resolution=None, y_resolution=None):
    """Creates flat file by Rasterizing layer.

    Args:
      input_shp: Input shapefile path.
      output_tif: Path to output raster tif.
      output_flat_shp: Output flat shapefile after polygonize the raster.
      width: Width of raster. Default value is 70,000. This is calculated
      to burn (Xmax - Xmin)/70000 this size pixel in width and
      (Ymax -Ymin)/35000 pixel size in height.
      height: Output raster height in pixel.
      attribute: Attribute value to be used in rasterization.
      x_resolution: Resolution of raster in degrees in x dimension.
      y_resolution: Resolution of raster in degrees in y dimension.

    Raises:
      Error: If unable to rasterize layer.
    """
    ds = ogr.Open(input_shp)
    layer = ds.GetLayer()
    xmin, xmax, ymin, ymax = layer.GetExtent()
    driver = gdal.GetDriverByName('GTiff')
    raster_dst = driver.Create(output_tif, width, height,
                               1, gdal.GDT_Int32, options=['COMPRESS=LZW'])
    raster_dst.GetRasterBand(1).Fill(0)
    raster_dst.GetRasterBand(1).SetNoDataValue(0)
    x_res = x_resolution or (xmax - xmin) / width
    y_res = y_resolution or (ymax - ymin) / height
    logging.info('Creating raster of x resolution: %s and y resolution: %s',
                 x_res, y_res)
    raster_transform = [xmin, x_res, 0, ymax, 0, -y_res]
    raster_dst.SetGeoTransform(raster_transform)
    srs_info = osr.SpatialReference()
    srs_info.SetWellKnownGeogCS('WGS84')
    raster_dst.SetProjection(srs_info.ExportToWkt())
    logging.info('Rasterizing layer....')
    # Blendrank value used as pixel value in rasterization. This value is
    # is really high which fits under 32bit.
    result = gdal.RasterizeLayer(
        raster_dst, [1], layer,
        burn_values=[_MAX_VALUE_32BIT],
        options=['ALL_TOUCHED=TRUE', 'ATTRIBUTE=%s' % attribute])
    if result:
      raise Error('Rasterization failed for layer: %s' % input_shp)

    src_band = raster_dst.GetRasterBand(1)
    logging.info('Creating data source ...')
    dst_ds = self._driver.CreateDataSource(output_flat_shp)
    dst_layer1 = dst_ds.CreateLayer(
        'flat_file_intermediate',
        geom_type=layer.GetLayerDefn().GetGeomType(),
        srs=layer.GetSpatialRef())
    fd = ogr.FieldDefn(attribute, ogr.OFTInteger)
    fd.SetWidth(15)
    dst_layer1.CreateField(fd)
    logging.info('Polygonizing layer....')

    gdal.Polygonize(src_band, src_band.GetMaskBand(), dst_layer1, 0,
                    ['8CONNECTED=8'])

  def GetFlatFileWithAllAttribute(self, raw_shp, flat_intmd_shp,
                                  shp_data_dict=None,
                                  burn_attribute='unique_id'):
    """Adds attribute to flat file.

    Args:
      raw_shp: Shapefile to read attributes.
      flat_intmd_shp: Path Intermediate shapefile.
      shp_data_dict: A dictionary mapping blendrank to attributes.
      If 'None' read attribute from raw shapefile.
      burn_attribute: Attribute field to be used for rasterize value.

    Raises:
      Error: If empty geometries are found in intermediate flat file.
    """
    dst_raw = ogr.Open(raw_shp)
    dst_layer_raw = dst_raw.GetLayer()

    dst_intmd_flat = ogr.Open(flat_intmd_shp)
    dst_layer_intmd = dst_intmd_flat.GetLayer()

    dst_final_flat = self._driver.CreateDataSource(self._output_flat_file)
    dst_layer_flat = dst_final_flat.CreateLayer(
        'flat_file_shp',
        geom_type=dst_layer_intmd.GetLayerDefn().GetGeomType(),
        srs=dst_layer_intmd.GetSpatialRef())
    attribute_list = []
    for i in range(dst_layer_raw.GetLayerDefn().GetFieldCount()):
      dst_layer_flat.CreateField(dst_layer_raw.GetLayerDefn().GetFieldDefn(i))
      attribute_list.append(
          dst_layer_raw.GetLayerDefn().GetFieldDefn(i).GetName())

    attribute_dict = {}
    if not shp_data_dict:
      for feature in dst_layer_raw:
        blend_rank = feature.GetField(burn_attribute)
        value_dict = {}
        for attribute in attribute_list:
          value_dict[attribute] = feature.GetField(attribute)
        attribute_dict[blend_rank] = value_dict
      shp_data_dict = attribute_dict

    flat_geom_dict_intmd = {}
    for feature in dst_layer_intmd:
      blendrank = feature.GetFieldAsInteger(burn_attribute)
      geom = feature.GetGeometryRef()
      # Rasterization process could create self intersecting polygons.
      # Self intersecting polygons boundary crosses itself and creates
      # problem in further processing. To fix these self intersecting
      # polygons, call CloseRings on the polygon geometry object.
      geom.CloseRings()
      if blendrank in flat_geom_dict_intmd:
        flat_geom_dict_intmd[blendrank].append(geom.ExportToWkb())
      else:
        flat_geom_dict_intmd[blendrank] = [geom.ExportToWkb()]
    if not flat_geom_dict_intmd:
      raise Error('No geometry found in intermediate flat file: %s'
                  % flat_intmd_shp)

    flat_feature = ogr.Feature(dst_layer_flat.GetLayerDefn())
    for blendrank in flat_geom_dict_intmd:
      geom_list = flat_geom_dict_intmd[blendrank]
      if len(geom_list) == 1:
        geom_wkb = geom_list[0]
      elif len(geom_list) > 1:
        try:
          geom_wkb_list = [wkb.loads(geom_wkb) for geom_wkb in geom_list]
          geom_shapely = ops.cascaded_union(
              [geom if geom.is_valid else geom.buffer(0)
               for geom in geom_wkb_list])
          geom_wkb = geom_shapely.wkb
        except ValueError as err:
          print('Found error: %s for %s: %s',
                        err, burn_attribute, blendrank)
          continue
      else:
        print('Empty geometry found for %s: %s',
                      burn_attribute, blendrank)
        continue
      poly_area = CalculateArea(geom_wkb)
      geom = ogr.CreateGeometryFromWkb(geom_wkb)

      flat_feature.SetGeometry(geom)
      for attribute in attribute_list:
        try:
          if attribute == 'area':
            flat_feature.SetField(
                attribute, poly_area)
          else:
            flat_feature.SetField(
                attribute, shp_data_dict[blendrank][attribute])
        except KeyError as err:
          print('Got error: %s for unknown %s: %s',
                        attribute, blendrank, err)
      dst_layer_flat.CreateFeature(flat_feature)


def CalculateArea(polygon_wkb):
  """Calculate are in Sqkm for given input wkt.

  To calcualte area of polygon from WKT which has WGS84 and lat long
  spatail reference is reproject in Cylindrical Equal-Area (EPSG:3410)
  projection.

  Args:
    polygon_wkb: Polygon well known binary.

  Returns:
    area: Area of polygon in SQKM.
  """
  src_sr = osr.SpatialReference()
  src_sr.ImportFromEPSG(_SRS_EPSG)
  dest_sr = osr.SpatialReference()
  dest_sr.ImportFromProj4(_CYLINDRICAL_EQUAL_AREA_PROJ4_DEF)
  sr_trans = osr.CoordinateTransformation(src_sr, dest_sr)
  new_geom = ogr.CreateGeometryFromWkb(polygon_wkb)
  geom = new_geom
  geom.Transform(sr_trans)
  area = geom.Area() / (1000 * 1000)  # SQKM
  return area


def GetTempFilePath(file_name):
  """Returns a unique temp file path for the given file name.

  Args:
    file_name: str, name of the file.
  Returns:
    a temporary path.
  """
  tempdir = tempfile.mkdtemp()
  return os.path.join(tempdir, file_name)


def main():
  parser = ArgumentParser(description = 'create flat file')
  parser.add_argument('--input_shp', required=True, help = "Path to input SHP file, in WGS84.")
  parser.add_argument('--flat_output_shp', required=True, help = "Path to flattened output SHP file.")
  parser.add_argument('--x_resolution', type=float, default=.001,
                    help="Resolution of raster in degrees in x direction.")
  parser.add_argument('--y_resolution', type=float, default=.001,
                    help="Resolution of raster in degrees in y direction.")
  parser.add_argument('--unique_id_field', type=str, default='unique_id',
                    help="Unique id field in input SHP.")
  parser.set_defaults(feature=True)
  args = parser.parse_args()

  # set up temp files.
  flat_int_coverage_path = GetTempFilePath('flat_coverage_int.shp')
  flat_coverage_tif = os.path.join(os.path.dirname(flat_int_coverage_path),
                                   'flat_coverage.tif')

  ff = FlatFile(args.flat_output_shp)
  ff.FlatLayerFromRasterize(
      args.input_shp,
      flat_coverage_tif,
      flat_int_coverage_path,
      attribute=args.unique_id_field,
      x_resolution=args.x_resolution,
      y_resolution=args.y_resolution)
  logging.info('Created %s', flat_int_coverage_path)
  ff.GetFlatFileWithAllAttribute(args.input_shp,
                                 flat_int_coverage_path,
                                 burn_attribute=args.unique_id_field)
  logging.info('Created %s', args.flat_output_shp)

if __name__ == "__main__":
  main()
