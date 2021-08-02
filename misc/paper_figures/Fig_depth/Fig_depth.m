%
% Saved for reference, but use the python script - see Fig_depth.py
%

load Fig_depth.h5

% Basic
figure(1)
plot (Frame001.autoalign.box_depth.box73, Frame001.signal.postproc.boxes.box73, '.')
xlabel('Depth (mm)')
ylabel('ISH signal')

% Histogrammed
figure(2)
