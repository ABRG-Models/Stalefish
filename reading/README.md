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

I tend to use GNU Octave for a quick look at a file, and python for analysis.
