# Reading the HDF5 data

This directory contains some help on getting access to the data file
written out by the program.

The data is written in "Hierarchical Data Format, version 5", or just
"HDF5", which is a widely-used, standard format.

There are several ways to read HDF5 files and get access to the data
they store, and you can choose the one that suits you best. Here's an
incomplete list:

* Python - using the h5py module or pandas.read_hdf()
* MATLAB - using the hdf5read function
* GNU Octave - It's as easy as "load myhdf5.h5"
* HDFView: [Download HDFView](https://portal.hdfgroup.org/display/support/Download+HDFView)
* The h5dump command, which can lead to a lot of text output, but it's
  all human readable.
* C or C++ - HDF5 has a C API, and morphologica has a wrapper class
  around the C API (See morph::HdfData).

I tend to use GNU Octave for a quick look at a file, and python or
morph::HdfData for analysis. Look at the README.md files in the octave
and python subdirectories for examples.

## Data layout

While the HDF5 format is standard, the meaning of the data is specific to
the Stalefish application, so here's some reference information on
what you'll find in the HDF5 files.

### Root level

At the root level of the data file, you'll find lots of 'Frames'. A
frame is the information about a single brain slice preparation. It
includes the filename of the image that was used, the points that the
user set to define the curve, the parameters for the gene expression
sampling bins and so
on. So, at the root level there are several data structures called
'Frame001', 'Frame002' and so on, as well as an integer 'nframes'
which makes it easy to write a loop through all the frames in an
analysis script.

 * Frame001, Frame002 etc: One Frame object for each brain slice.
 * nframes: The number of Frame objects in this data file.

### Frame level

Each frame contains:

 * boxes_pixel0, boxs_pixel1, boxs_pixel2, etc, each of which is a
   list of the pixel values for the "sampling boxes" - the yellow
   boxes that are visible in the user interface main window (and which
   are black on the signal window). Range: 0 to 255.

 * boxes_signal0, boxs_signal1, boxs_signal2, etc, each of which is a
   list of the signal values for the sampling boxes. The signal is
   derived from the pixel values, with compensation for any systematic
   variation in the background illumination. Range: Approx 0 to 1
   (floating point).

 * FIXME: Need corresponding lists of coordinates in pixels and in
   autoaligned/lmaligned coordinates.

 * class - a structure containing all the data required to re-open a
   Stalefish project, some of which may be useful in analysis. See
   below for a description of the objects which you will find inside
   class.

 * fitted - a list of the fitted points for the curve fitted to the
   user-supplied points. What coordinate system? Pixels? FIXME

 * autoalign_computed A true/false flag. If true (>0) then
   the brain slice auto-alignment algorithm based on the curves was
   computed and stored in fitted_autoaligned.

 * fitted_autoaligned - a list of the fitted points for the curve fitted to the
   user-supplied points. Auto-aligned coordinates (mm)

 * lmalign_computed A true/false flag. If true (>0) then
   the brain slice landmark alignment algorithm based on user-supplied
   landmarks was
   computed and stored in fitted_lmaligned

 * fitted_lmaligned - a list of the fitted points for the curve fitted to the
   user-supplied points. Landmark-aligned 2D coordinates (mm)

 * LM_autoaligned - 3d coordinates of the landmarks for this slice,
   transformed as per the other _autoaligned points.

 * LM_lmaligned - 3d coordinates of the landmarks for this slice,
   transformed as per the other landmark aligned points.

 * box_pixel_means - The mean pixel value for each box. Can be re-created by
   computing the means of boxes_pixel0, boxes_pixel1, etc.

 * box_signal_means - The mean signal

 * box_signal_means_autoscaled - ??

 * sbox_centers_autoaligned - The "surface box" centers, as 3D coordinates. A
   'surface box' is the box which is a long as the corresponding bin
   is wide, and as wide as the slice is thick. This used the
   2D coordinates from fitted_autoalign to construct these 3D
   coordinates.

 * sboxes_autoaligned - The surface boxes, all four corners for the
   auto-aligned slices. Each element of this list
   contains 12 numbers, which specify 4 coordinates in 3D space; the
   corners of the surface boxes.

* sbox_centers_lmaligned - The "surface box" centers, as 3D coordinates. A
   'surface box' is the box which is a long as the corresponding bin
   is wide, and as wide as the slice is thick. This used the
   2D coordinates from fitted_autoalign to construct these 3D
   coordinates.

 * sboxes_lmaligned - The surface boxes, all four corners for the
   landmark-aligned slices.

 * sbox_centers_scaled - The surface box centers scaled from pixels to
   mm, but with no other alignments. These are saved for debugging
   purposes.

 * sboxes_scaled - The four corners of the
   scaled-but-otherwise-untranslated surface boxes.

 * sbox_linear_distance - The linear distance along the curve of each
   **autoaligned** surface box. These distances are centered so that the
   middle surface box should have a value of 0 for the
   sbox_linear_distance.

 * nboxes: The number of box objects.

 * freehand0, freehand1, etc, each of which is a list of the pixel
   values in freehand drawn regions (drawn in 'freehand' mode)

 * freehand_means The mean pixel value in each freehand drawn region.

 * FIXME: Need to store the freehand pixel coordinates.

 * nfreehand The number of freehand drawn regions.

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
 * FLE000, FLE001, etc: 2D coordinates of the pixels enclosed in each
  freehand loop.
 * FLE_n: The number of freehand loops
