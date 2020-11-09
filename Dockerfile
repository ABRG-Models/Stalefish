FROM gcc:latest

RUN set -ex;                                                                      \
    apt-get update;                                                               \
    apt-get install -y cmake libzmq3-dev; 										  \
    apt-get install cmake gcc g++;												  
    
ENV CC /usr/bin/clang
ENV CXX /usr/bin/g++
ENV OPENCV_VERSION='3.4.2' DEBIAN_FRONTEND=noninteractive

RUN cd /opt && \
  wget https://github.com/opencv/opencv/archive/${OPENCV_VERSION}.zip && \
  unzip ${OPENCV_VERSION}.zip && \
  rm -rf ${OPENCV_VERSION}.zip

RUN mkdir -p /opt/opencv-${OPENCV_VERSION}/build && \
  cd /opt/opencv-${OPENCV_VERSION}/build && \
  cmake \
    -D BUILD_DOCS=OFF \
    -D BUILD_EXAMPLES=OFF \
    -D BUILD_opencv_apps=OFF \
    -D BUILD_opencv_python2=OFF \
    -D BUILD_opencv_python3=OFF \
    -D BUILD_PERF_TESTS=OFF \
    -D BUILD_SHARED_LIBS=OFF \ 
    -D BUILD_TESTS=OFF \
    -D CMAKE_BUILD_TYPE=RELEASE \
    -D ENABLE_PRECOMPILED_HEADERS=OFF \
    -D FORCE_VTK=OFF \
    -D WITH_FFMPEG=OFF \
    -D WITH_GDAL=OFF \ 
    -D WITH_IPP=OFF \
    -D WITH_OPENEXR=OFF \
    -D WITH_OPENGL=OFF \ 
    -D WITH_QT=OFF \
    -D WITH_TBB=OFF \ 
    -D WITH_XINE=OFF \ 
    -D BUILD_JPEG=ON  \
    -D BUILD_TIFF=ON \
    -D BUILD_PNG=ON \
  .. && \
  make -j$(nproc) && \
  make install && \
  rm -rf /opt/opencv-${OPENCV_VERSION}

COPY . /usr/src/test

WORKDIR /usr/src/test

RUN mkdir build

RUN cd build

RUN cmake .

RUN make

CMD ["./build/src/stalefish ./data/testimg.json"]