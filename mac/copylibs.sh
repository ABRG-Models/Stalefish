#/bin/bash

# cmake puts rpaths into the executables by default
RPATHS_IN_EXE=1

# Create Stalefish.app tree
mkdir -p Stalefish.app/Contents/MacOS
mkdir -p Stalefish.app/Contents/Resources

# Copy in the icons
cp sf_icons.icns Stalefish.app/Contents/Resources/

# Copy the freshly built executable(s)
PROGRAM_LIST="stalefish sfview sfresave sfgetjson"
for j in ${PROGRAM_LIST}; do
    cp ../build/src/${j} Stalefish.app/Contents/MacOS/
done

STALEFISH_LIBS="libhdf5.103.dylib libarmadillo.9.dylib libjsoncpp.24.dylib libopencv_dnn.4.5.dylib libopencv_gapi.4.5.dylib libopencv_highgui.4.5.dylib libopencv_ml.4.5.dylib libopencv_objdetect.4.5.dylib libopencv_photo.4.5.dylib libopencv_stitching.4.5.dylib libopencv_video.4.5.dylib libopencv_videoio.4.5.dylib libopencv_imgcodecs.4.5.dylib libopencv_calib3d.4.5.dylib libopencv_features2d.4.5.dylib libopencv_flann.4.5.dylib libopencv_imgproc.4.5.dylib libopencv_core.4.5.dylib"

for i in ${STALEFISH_LIBS}; do

    # This copies the versions on Seb's Mac into the package
    cp /usr/local/lib/${i} Stalefish.app/Contents/MacOS/

    # This would change the links to the libraries in the executable, if,
    # in the executable rpaths weren't already used, but they ARE in the
    # STALEFISH_LIBS
    if [ ${RPATHS_IN_EXE} -ne 1 ]; then
        for j in ${PROGRAM_LIST}; do
            install_name_tool -change /usr/local/lib/${i} @rpath/${i} Stalefish.app/Contents/MacOS/${j}
        done
    fi
done

# Extra libs required by sfview
cp /usr/local/lib/libomp.dylib Stalefish.app/Contents/MacOS/

cp /usr/local/lib/libpopt.0.dylib Stalefish.app/Contents/MacOS/
# My build process doesn't make sfview access libpopt with RPATH, so fix that
install_name_tool -change /usr/local/lib/libpopt.0.dylib @rpath/libpopt.0.dylib Stalefish.app/Contents/MacOS/sfview

cp /opt/X11/lib/libfreetype.6.dylib Stalefish.app/Contents/MacOS/
# My build process doesn't make sfview access ft with RPATH, so fix that too
install_name_tool -change /opt/X11/lib/libfreetype.6.dylib @rpath/libfreetype.6.dylib Stalefish.app/Contents/MacOS/sfview

# Ensure the RPATH is set for our bundled libraries
for j in ${PROGRAM_LIST}; do
    # @loader_path is the path to the folder containing the program (${j})
    install_name_tool -add_rpath @loader_path Stalefish.app/Contents/MacOS/${j}
done

# Create the Info.plist file
cat > Stalefish.app/Contents/Info.plist <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
	<key>CFBundleExecutable</key>
	<string>stalefish</string>
	<key>CFBundleIconFile</key>
	<string>sf_icons.icns</string>
	<key>CFBundleIdentifier</key>
	<string>com.academic.sebjames</string>
	<key>NSHighResolutionCapable</key>
	<true/>
</dict>
</plist>
EOF
