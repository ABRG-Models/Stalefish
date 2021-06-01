#/bin/bash

# Create Stalefish.app tree
mkdir -p Stalefish.app/Contents/MacOS
mkdir -p Stalefish.app/Contents/Resources

# Copy in the icons
cp sf_icons.icns Stalefish.app/Contents/Resources/

# Copy the freshly built executable(s)
PROGRAM_LIST="stalefish"
for j in ${PROGRAM_LIST}; do
    cp ../build/src/${j} Stalefish.app/Contents/MacOS/
done

for i in libhdf5.103.dylib libarmadillo.9.dylib libjsoncpp.24.dylib libopencv_dnn.4.5.dylib libopencv_gapi.4.5.dylib libopencv_highgui.4.5.dylib libopencv_ml.4.5.dylib libopencv_objdetect.4.5.dylib libopencv_photo.4.5.dylib libopencv_stitching.4.5.dylib libopencv_video.4.5.dylib libopencv_videoio.4.5.dylib libopencv_imgcodecs.4.5.dylib libopencv_calib3d.4.5.dylib libopencv_features2d.4.5.dylib libopencv_flann.4.5.dylib libopencv_imgproc.4.5.dylib libopencv_core.4.5.dylib; do

    # This copies the versions on Seb's Mac into the package
    cp /usr/local/lib/${i} Stalefish.app/Contents/MacOS/

    # This changes the executable
    for j in ${PROGRAM_LIST}; do
        install_name_tool -change /usr/local/lib/${i} @executable_path/${i} Stalefish.app/Contents/MacOS/${j}
    done
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
