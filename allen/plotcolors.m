% The expressing group
b = unique(a, "rows");
% The non-expressing group
o = unique(n, "rows");

len2 = 100

lenb = size(b)(1)
rndIDX = randperm(lenb);
% ex is expressing
ex = b(rndIDX(1:len2), :);

leno = size(o)(1)
rndIDX = randperm(leno);
% ne is non-expressing
ne = o(rndIDX(1:len2), :);

% Note: marker face colour is the colour the pixel has on the ISH
% image.
figure (1)
clf
hold on
ex_col = [111,0,207]./255.0;
nex_col = [12,140,84]./255.0;
for ii = 1:len2
    ex_x = [ex(ii,1)];%b
    ex_y = [ex(ii,2)];%g
    ex_z = [ex(ii,3)];%r
    plot3(ex_x,ex_y,ex_z,'o','markeredgecolor',ex_col,'markersize', 10, 'markerfacecolor', flip(ex(ii,:)./255.0));
    ne_x = [ne(ii,1)];
    ne_y = [ne(ii,2)];
    ne_z = [ne(ii,3)];
    plot3(ne_x,ne_y,ne_z,'o','markeredgecolor',nex_col,'markersize', 10, 'markerfacecolor', flip(ne(ii,:)./255.0));
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
ex_col = [111,0,207]./255.0;
nex_col = [12,140,84]./255.0;
blue = [];
greenred = [];
for ii = 1:len2
    ex_x = [ex(ii,1)];%b
    ex_y = [ex(ii,2)];%g
    ex_z = [ex(ii,3)];%r
    plot3(ex_x,ex_y,ex_z,'o','markeredgecolor',ex_col,'markersize', 10, 'markerfacecolor', flip(ex(ii,:)./255.0));
    blue = [blue; ex_x, 1];
    greenred = [greenred; ex_y, ex_z];
end

plot3([0 255], [0 0], [0 0], 'b-');
plot3([0 0], [0 255], [0 0], 'g-');
plot3([0 0], [0 0], [0 255], 'r-');

xlabel ('Blue channel')
ylabel ('Green channel')
zlabel ('Red channel')

% Linear regression. Regress red/green against blue
thefit = blue\greenred;
greenred_calc = (thefit(1,:) .* blue(:,1)) + thefit(2,:);
% Plot best-fit line
plot3 (blue(:,1), greenred_calc(:,1), greenred_calc(:,2), 'k-')

% Add Lydia Ng's luminosity line:
t = [0:0.01:1];
rt = 0.21 .* t .*255;
gt = 0.72 .* t .*255;
bt = 0.07 .* t .*255;
ng_col = [ 0.21, 0.72, 0.07 ];
plot3(bt, gt, rt, '.-', 'markerfacecolor', ng_col, 'color', ng_col);

% Add the fit line to fig 1 too
figure(1)
plot3 (blue(:,1), greenred_calc(:,1), greenred_calc(:,2), 'k.-');

% Add Lydia Ng's luminosity line:
plot3(bt, gt, rt, 'v-', 'markerfacecolor', ng_col, 'color', ng_col);

%
% Translating and rotating points so they lie on the x axis
%
figure (3);
clf;
hold on;

% Blue input, evenly spaced
b_in = [0:1:255]';

% Save the translation offset from the fit
trans_offset = thefit(2,:);
trans_offset3d = [0.0, thefit(2,:)];

% Resulting calculated values based on fit
rg_calc = (thefit(1,:) .* b_in) + trans_offset;
% Translation - just subtract the offset in the fit.
rg_trans = rg_calc - trans_offset;
% Plot the translated fit
plot3(b_in, rg_trans(:,1), rg_trans(:,2), 'color', 'm');

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
Ay = [cos(theta), 0, sin(theta); ...
      0, 1, 0; ...
      -sin(theta), 0, cos(theta)];

Az = [cos(-phi), -sin(-phi), 0; ...
      sin(-phi), cos(-phi), 0; ...
      0, 0, 1];

% Create original points matrix:
original_points = [blue(:,1), greenred];
% Create original points with the translation, too:
gr_trans = greenred - trans_offset;
translated_points = [blue(:,1), gr_trans];

% Create a single rotation matrix A from Ay and Az
A = Ay * Az;
% These translated, rotated points have 2 rotations applied:
trans_rot_points = (A*translated_points')';
plot3 (trans_rot_points(:,1),trans_rot_points(:,2),trans_rot_points(:,3),'ok');

plot3([0 255], [0 0], [0 0], 'b-');
plot3([0 0], [0 255], [0 0], 'g-');
plot3([0 0], [0 0], [0 255], 'r-');

xlabel ('Blue channel')
ylabel ('Green channel')
zlabel ('Red channel')
title ('Transformed points')

%
% Next job: find first two principle components in the greenred plane and use
% these to constrain an ellipse that will select those pixels whose
% colour is close enough to be considered to be expressing.
%
figure(4)
clf
X = trans_rot_points(:,2:3);
plot (X(:,1),X(:,2),'o');

% PCA.
% My data are already centered around 0, so no need to subtract the mean:
mu = mean(X);
%Xm = bsxfun(@minus, X, mu)
C = cov(X);
[V,D] = eig(C);

% sort eigenvectors desc
[D, i] = sort(diag(D), 'descend');
V = V(:,i);
% V are eigenvectors, D eigen values - they're like 'variances'
hold on
scale = 2.*sqrt(D(1));
pc1_x = mu(1) + scale * V(1,1);
pc1_y = mu(2) + scale * V(2,1);
pc1 = line([0, pc1_x], [0, pc1_y]);
scale2 = 2.*sqrt(D(2));
pc2_x = mu(1) + scale2 * V(1,2);
pc2_y = mu(2) + scale2 * V(2,2);
pc2 = line([0, pc2_x], [0, pc2_y]);
set(pc1, 'color', [1 0 0], "linestyle", "--");
set(pc2, 'color', [0 1 0], "linestyle", "--");
% Draw an ellipse.
gamma_rad = atan2(pc1_y, pc1_x);
gamma = 360 * gamma_rad / (2*pi);
% This from Octave's geometry package
drawEllipse(0, 0, 2.*sqrt(D(1)), 2.*sqrt(D(2)), gamma);
% The key is that the principle axis lengths of the ellipse

% Un-rotate the ellipse:
drawEllipse(0, 0, 2.*sqrt(D(1)), 2.*sqrt(D(2)), 0);
% Rotate the points, too
A_2d = [cos(-gamma_rad), -sin(-gamma_rad); sin(-gamma_rad), cos(-gamma_rad)];
r_2d = (A_2d * X')';
plot (r_2d(:,1), r_2d(:,2),'om');

% I've now got a 3rd rotational transformation, so Let's apply this
Ax = [1, 0, 0;...
      0, cos(-gamma_rad), -sin(-gamma_rad); ...
      0, sin(-gamma_rad), cos(-gamma_rad)];

figure(5)
clf
hold on;

% One way to compute the fully transformed/rotated points:
%fully_transformed = (Ax*trans_rot_points')';
%plot3 (fully_transformed(:,1),fully_transformed(:,2),fully_transformed(:,3),'ok');

% Note that the output of this script can now be a single rotation
% matrix (A) along with the two axes of the ellipse to form the
% equation to categorize a color as expressing or not expressing.
A = Ax * Ay * Az;
fully_transformed = (A*translated_points')';
plot3 (fully_transformed(:,1), fully_transformed(:,2), fully_transformed(:,3), '.k');

% Let's make a one-shot transformation matrix
oneshot_figured = 0
if oneshot_figured
    % 4 vector version of orig points:
    original_points4vec = [original_points, ones(length(original_points),1)];

    % Trans matrix
    T = [1 0 0; 0 1 0; 0 0 1];
    T = [T, trans_offset3d'];
    T = [T; 0 0 0 1]

    translated_points2 = (T * original_points4vec')';

    % Rot matrxi
    R = [A, [0;0;0]];
    R = [R; 0 0 0 1]

    fully_transformed2 = (R*translated_points2')';

    % Together:
    AT = [A, -trans_offset3d'];
    AT = [AT; 0 0 0 1]
    fully_transformed3 = (AT*original_points4vec')';

    plot3 (fully_transformed2(:,1),fully_transformed2(:,2),fully_transformed2(:,3),'ob');
    plot3 (fully_transformed3(:,1),fully_transformed3(:,2),fully_transformed3(:,3),'og');
end

% For comparison:
%plot3 (trans_rot_points(:,1),trans_rot_points(:,2),trans_rot_points(:,3),'.m');

% Let's plot the transformed non-expressing group in red:
n_expr = (A * (ne-trans_offset3d)')';
plot3 (n_expr(:,1),n_expr(:,2),n_expr(:,3),'*r');

%% Lets draw the ellipse on Fig 5.
ellip_major = 2.*sqrt(D(1));
ellip_minor = 2.*sqrt(D(2));
% y = +/- ellip_minor .* sqrt (1 - x*x/(ellip_major*ellip_major))
% What's the range of X? It's +- ellip_major
ellip_x = [-ellip_major:1:ellip_major]';
ellip_y_plus = + ellip_minor .* sqrt (1 - ellip_x.*ellip_x/(ellip_major*ellip_major));
ellip_y_minus = - ellip_minor .* sqrt (1 - ellip_x.*ellip_x/(ellip_major*ellip_major));
ellip_plus = [300.*ones(size(ellip_x),1), ellip_x, ellip_y_plus];
ellip_minus = [300.*ones(size(ellip_x),1), ellip_x, ellip_y_minus];
plot3(ellip_plus(:,1), ellip_plus(:,2), ellip_plus(:,3), 'ob-')
plot3(ellip_minus(:,1), ellip_minus(:,2), ellip_minus(:,3), 'ob-')
%% Done drawing ellipse

%% Axes
plot3([0 255], [0 0], [0 0], 'b-');
plot3([0 0], [0 255], [0 0], 'g-');
plot3([0 0], [0 0], [0 255], 'r-');
xlabel ('Blue channel')
ylabel ('Green channel')
zlabel ('Red channel')
title ('Final transformed points')

% Last question. Is a point inside the ellipse? Plot, and colour
% accordingly
%
% For a given x,y(red,green) in transformed colour space, do the
% sum:
inout = []
values = []
for f = fully_transformed'
    erad = ((f(2).*f(2))./(ellip_major.*ellip_major)) + ((f(3).*f(3))./(ellip_minor.*ellip_minor));
    if (erad > 1)
        %disp('erad > 1')
        % Outside, so value is 0
        inout = [inout 0];
        values = [values 0];
    else
        %disp('erad <= 1')
        inout = [inout 1];
        values = [values f(1)];
    end
end

% Apply linear transformation to values:
luminosity_cutoff = 250;
luminosity_factor = -1;
tvalues = luminosity_cutoff + (values(inout==1) .* luminosity_factor);
tvalues(tvalues<0) = 0;
plot3(fully_transformed(inout==1,1), zeros(length(values(inout==1)),1), tvalues, 'go');

% Process each number through the eqn of the ellipse to get a vlaue
% < or > 1.

printf ('================ RESULTS ==================\n')
printf ('3D translation:\n')
trans_offset3d
printf ('3D rotation matrix:\n')
A
printf ('major & minor axes of ellipse in the (transformed )red-green plane:\n')
ellip_major
ellip_minor
printf ('Luminosity factor and cutoff\n')
luminosity_factor
luminosity_cutoff
printf ('===========================================\n')
