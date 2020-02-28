b = unique(a, "rows");
o = unique(n, "rows");

len2 = 200

lenb = size(b)(1)
rndIDX = randperm(lenb);
c = b(rndIDX(1:len2), :);

leno = size(o)(1)
rndIDX = randperm(leno);
p = o(rndIDX(1:len2), :);

% Note: marker face colour is the colour the pixel has on the ISH
% image.
figure (1)
clf
hold on
ex_col = [111,0,207]./255.0
nex_col = [12,140,84]./255.0
for ii = 1:len2
    ex_x = [c(ii,1)];%b
    ex_y = [c(ii,2)];%g
    ex_z = [c(ii,3)];%r
    plot3(ex_x,ex_y,ex_z,'o','markeredgecolor',ex_col,'markersize', 10, 'markerfacecolor', flip(c(ii,:)./255.0));
    ne_x = [p(ii,1)];
    ne_y = [p(ii,2)];
    ne_z = [p(ii,3)];
    plot3(ne_x,ne_y,ne_z,'o','markeredgecolor',nex_col,'markersize', 10, 'markerfacecolor', flip(p(ii,:)./255.0));
end

plot3([0 255], [0 0], [0 0], 'b-');
plot3([0 0], [0 255], [0 0], 'g-');
plot3([0 0], [0 0], [0 255], 'r-');

xlabel ('Blue channel')
ylabel ('Green channel')
zlabel ('Red channel')

figure (2)
clf
hold on
ex_col = [111,0,207]./255.0
nex_col = [12,140,84]./255.0
blue = [];
redgreen = [];
for ii = 1:len2
    ex_x = [c(ii,1)];%b
    ex_y = [c(ii,2)];%g
    ex_z = [c(ii,3)];%r
    plot3(ex_x,ex_y,ex_z,'o','markeredgecolor',ex_col,'markersize', 10, 'markerfacecolor', flip(c(ii,:)./255.0));
    blue = [blue; ex_x, 1];
    redgreen = [redgreen; ex_y, ex_z];
end

plot3([0 255], [0 0], [0 0], 'b-');
plot3([0 0], [0 255], [0 0], 'g-');
plot3([0 0], [0 0], [0 255], 'r-');

xlabel ('Blue channel')
ylabel ('Green channel')
zlabel ('Red channel')

% Linear regression. Regress red/green against blue
thefit = blue\redgreen;
redgreen_calc = (thefit(1,:) .* blue(:,1)) + thefit(2,:);
% Plot best-fit line
plot3 (blue(:,1), redgreen_calc(:,1), redgreen_calc(:,2), 'k-')

% Add Lydia Ng's luminosity line:
t = [0:0.01:1];
rt = 0.21 .* t .*255;
gt = 0.72 .* t .*255;
bt = 0.07 .* t .*255;
ng_col = [ 0.21, 0.72, 0.07 ];
plot3(bt, gt, rt, '.-', 'markerfacecolor', ng_col, 'color', ng_col);

% Add the fit line to fig 1 too
figure(1)
plot3 (blue(:,1), redgreen_calc(:,1), redgreen_calc(:,2), 'k-')

% Add Lydia Ng's luminosity line:
plot3(bt, gt, rt, 'markerfacecolor', ng_col, 'color', ng_col);

%
% Translating and rotating points so they lie on the x axis
%

% Translation - just subtract the offset in the fit.
redgreen_trans = redgreen_calc - thefit(2,:);
plot3(blue(:,1), redgreen_trans(:,1), redgreen_trans(:,2), 'color', 'r');

% Rotation given by thefit(1,:)
x = 1;
yz = thefit(1,:) .* x;
y = thefit(1,1) .* x;
z = thefit(1,2) .* x;

% theta is rotation down to the x-y plane; phi is rotation about
% the z axis.
phi = atan2 (y, x);
theta = atan2 (z, sqrt(x*x + y*y));

% Create rotation matrices
Az = [cos(phi), -sin(phi), 0; ...
      sin(phi), cos(phi), 0; ...
      0,0,1];

Ay = [cos(theta), 0, sin(theta); ...
      0,1,0; ...
      -sin(theta) 0 cos(theta)];

A = Az*Ay;

Ay_back = [cos(-theta), 0, sin(-theta);...
           0,1,0;...
           -sin(-theta), 0, cos(-theta)];

Az_back = [cos(-phi),-sin(-phi),0;...
           sin(-phi),cos(-phi),0;...
           0,0,1];
A_back = Ay_back * Az_back;

points = [blue(:,1), redgreen];

tpoints = (Ay_back*points')';

figure(3)
hold on;
for ii = 1:len2
    ex_x = [tpoints(ii,1)];%b
    ex_y = [tpoints(ii,2)];%g
    ex_z = [tpoints(ii,3)];%r
    plot3(ex_x,ex_y,ex_z,'o','markeredgecolor',ex_col,'markersize', 10, 'markerfacecolor', [1,0,0]);
end

tpoints = (Az_back*points')';

figure(3)
hold on;
for ii = 1:len2
    ex_x = [tpoints(ii,1)];%b
    ex_y = [tpoints(ii,2)];%g
    ex_z = [tpoints(ii,3)];%r
    plot3(ex_x,ex_y,ex_z,'o','markeredgecolor',ex_col,'markersize', 10, 'markerfacecolor', [0,0,1]);
end

plot3([0 255], [0 0], [0 0], 'b-');
plot3([0 0], [0 255], [0 0], 'g-');
plot3([0 0], [0 0], [0 255], 'r-');

xlabel ('Blue channel')
ylabel ('Green channel')
zlabel ('Red channel')
