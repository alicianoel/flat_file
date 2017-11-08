# create_flat_file
Create flattened shapefile of non-overlapping polygons

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
  /path/to/flat_output_shp.shp 
  
  
