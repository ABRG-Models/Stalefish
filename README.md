# Stalefish

Stalefish is a tool for the analysis of 2D brain slice images, allowing the neurobiologist to reconstruct three dimensional gene expression surfaces from ISH slices. The 3D surfaces can then be digitally unwrapped to make 2D expression maps. A paper describing the process is under submission for publication.

Stalefish is written in C++ and will compile on Mac OS or Linux. Released versions are available as .dmg images or on the snapcraft.io store.

## Build dependencies

### Ubuntu 20.04+ Quick-start

Ubuntu 20.04 has apt-installable packages for all of Stalefish's dependencies! This should be a complete recipe:

```
sudo apt install build-essential cmake git wget  \
                 freeglut3-dev libglu1-mesa-dev libxmu-dev libxi-dev liblapack-dev \
                 libopencv-dev libarmadillo-dev libjsoncpp-dev libglfw3-dev \
                 libhdf5-dev libfreetype-dev libpopt-dev
```

### Dependencies on Mac OS or older Linux systems

#### morphologica dependencies

This program compiles with
[morphologica](https://github.com/ABRG-Models/morphologica), which is included as a git submodule. This means that Stalefish
needs to link to the morphologica-associated libraries armadillo, opencv, glfw, hdf5 and freetype. So head over
to the [morphologica Linux installation readme](https://github.com/ABRG-Models/morphologica/blob/main/README.install.linux.md) or the [Mac installation readme](https://github.com/ABRG-Models/morphologica/blob/main/README.install.mac.md) and follow the instructions to install the *dependencies* on your OS (you *don't* need to *build* morphologica).

#### libpopt

Stalefish also uses libpopt in the sfview tool, so build libpopt-1.18:
```
CPPFLAGS="-mmacosx-version-min=10.14" LDFLAGS="-mmacosx-version-min=10.14" ./configure --prefix=/usr/local
make
sudo make install
```
(Note I used -mmacosx-version-min=10.14 to build for a minimum Mac OS of Mojave. You don't need to do that; it's completely optional on a Mac and definitly leave it out on Linux!)

#### libomp on Mac

I used some OpenMP pragmas in sfview, so on a Mac, you need to make sure you have libomp. This is compiled from the llvm compiler source code.

```
git clone https://github.com/llvm/llvm-project.git
cd llvm-project/openmp
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/usr/local -DCMAKE_OSX_DEPLOYMENT_TARGET=10.14 ..
make omp
```
(An equivalent installation is not necessary on Linux)

## Compile Stalefish using cmake

```bash
# clone Stalefish and change into Stalefish dir:
git clone https://github.com/ABRG-Models/Stalefish.git
cd Stalefish
# Inside Stalefish dir, create a build dir:
mkdir build
# Change directory to build and do a standard cmake build
cd build
cmake ..
make
```

## Updating Stalefish

```
cd /path/to/Stalefish
git pull -p
cd build
make
```

## Running the program

run using e.g.,:
```
./build/src/stalefish ./data/testimg.json
```

The JSON file contains a list of the images to fit curves to, and some
other parameters.

Left-click points to which a curve should be fit. Bezier curves fit
best to 4 or 5 points; pressing space will create a new curve to fit;
the program will ensure they join nicely. Sliders give control over
the size and shape of the bins.

Press "h" to see help text.

## Analysing the data from the program

See the reading/ subdirectory and its README.md file.

## JSON parameters

Stalefish projects are created with a hand-written .json file which lists the images in your brain slice set, along with some additional information, such as the position of each slice in the stack, the slice thickness and so on. The parameters are listed here. Example json files can be found in the data/ directory - see testimg.json and vole_65_7E_id2_L23.json for templates.

* **pixels_per_mm** Set to the number of pixels per mm in the original image files.
* **thickness** The thickness of a brain slice (in mm), assuming they'll all be the same.
* **bg_blur_screen_proportion**: Float. Typically 1/6 (0.1667). The sigma for the Gaussian used to blur the image to get the overall background luminance is the framewidth in pixels multiplied by this number.
* **bg_blur_subtraction_offset**: Float. Range 0.0 to 255.0. Typically 180. A subtraction offset used when subtracting blurred background signal from image
* **colourmodel**: Enumerated. Options are *monochrome* (default) or *Allen*
* **colour_rot** Array related to Allen colour images
* **colour_trans** Array related to Allen colour images
* **ellip_axes** Array (of 2 numbers) related to Allen colour images
* **luminosity_cutoff**: Float. No longer used in monochrome mode.
* **luminosity_factor**: Float. No longer used in monochrome mode.
* **save_per_pixel_data**: Boolean. If true, save coordinates and signal value for every pixel in every sample box. Inflates .h5 file somewhat
* **save_auto_align_data**: Boolean. If true, save the slice auto-alignment location data
* **save_landmark_align_data**: Boolean. If true, save the slice landmark-alignment location data
* **save_frame_image**: Boolean. If true, save all images in .h5 file.
* **scaleFactor**: Float. By what factor should the images given in **slices** be scaled before they are shown in the UI. This is intended to allow you to scale down your images so that they fit within the resolution of your screen, or so that the data files generated by the program are kept in check. By scaling down the image size, the number of data points saved by the program is reduced, especially when **save_per_pixel_data** is set true
* **rotate_landmark_one**: Boolean. If true, and there are >1 landmark per slice, apply the *rotate slices about landmark 1* alignment procedure anyway. This rotational alignment is applied by default if there is ONLY 1 landmark per slice.
* **rotate_align_landmarks**: Boolean. If true, in "rotate about landmark 1 mode" align the other landmarks, instead of the curves.
* **slices**: Array of JSON objects specifying slice images **filename** and the slice's **x** position.
