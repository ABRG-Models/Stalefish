# Reading the HDF5 data

This directory contains some help on getting access to the data file
written out by the program.

The data is written in "Hierarchical Data Format, version 5", or just
HDF5, which is a global standard format.

There are several ways to read HDF5 files and get access to the data
they store, and you can choose the one that suits you best. Here's an
incomplete list:

* Python - using the h5py module or pandas.read_hdf()
* MATLAB - using the hdf5read function
* GNU Octave - It's as easy as "load myhdf5.h5"
* HDFView: [Download HDFView](https://portal.hdfgroup.org/display/support/Download+HDFView)

I tend to use GNU Octave for a quick look at a file, and python for
analysis. Look at the README.md files in the octave and python
subdirectories for examples.

## Data layout

So the format is standard, but the meaning of the data is specific to
the Stalefish application, so here's some reference information on
what you'll find in the HDF5 files.

### Root level

At the root level of the data file, you'll find lots of 'Frames'. A
frame is the information about a single brain slice preparation. It
includes the filename of the image that was used, the points that the
user set to define the curve, the parameters for the gene expression
sampling bins and so
on. So, at the root level there are several data structures called
'Frame000', 'Frame001' ann so on, as well as an integer 'nframes'
which makes it easy to write a loop through all the frames in an
analysis script.

 * Frame000, Frame001 etc: One Frame object for each brain slice.
 * nframes: The number of Frame objects in this data file.

### Frame level

Each frame contains:

 * box0, box1, box2, etc, each of which is a list of the pixel values
   for the bins.
 * class - a structure containing all the data required to re-open a
   Stalefish project, some of which may be useful in analysis
 * fitted - a list of the fitted points for the curve fitted to the
   user-supplied points.
 * fitted_offset - a list of the fitted points, offset by the relative
   centroids of the previous frame.
 * fitted_rotated - a list of the offset fitted points, but rotated by
   an amount that minimises the distances between the points on this
   curve and the points on the curve from the previous frame.
 * means - The mean pixel value for each box. Can be re-created by
   computing the means of box0, box1, etc.
 * means_autoscaled - The autoscaled version of mean (scaled from 0 to 1).
 * sbox_centers - The "surface box" centers, as 3D coordinates. A
   'surface box' is the box which is a long as the corresponding bin
   is wide, and as wide as the slice is thick.
 * sbox_linear_distance - The linear distance along the curve of each
   surface box. These distances are centered so that the middle
   surface box should have a value of 0 for the sbox_linear_distance.
 * sboxes - The surface box coordinates. Each element of this list
   contains 12 numbers, which specify 4 coordinates in 3D space; the
   corners of the surface boxes.
 * nboxes: The number of box objects.

### Class level

Each frame's class object contains:

 * P: Coordinates of the user points which have not yet be 'finalised'
  as a complete 'curve set'. These are the points which appear green
  in the application.
 * PP000, PP001, etc: 2D Coordinate 'curve sets' - the red and blue
  'user points'.
 * PP_n: The number of PP* curve sets.
 * binA: The distance from the fitted curve, along a normal to the curve, at which the bins start.
 * binB: The distance from the fitted curve, along a normal to the curve, at which the bins stop.
 * filename: The name of the image file of the brain slice
 * flags: Some flags which specify how the application should
  visualise the curve
 * idx: The frame's index number
 * layer_x: The position of the brain slice along the x axis.
 * nBinsTarg: The number of bins along the fitted curve.
 * pixels_per_mm: The scale for the brain slice image.
 * polyOrder: The order of the polynomial curve, if used.
 * pp_idx: A state variable for the program (not useful for analysis).
 * thickness: The thickness of the brain slice, as copied from the
  JSON config file.
