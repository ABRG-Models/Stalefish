b = unique(a, "rows");
o = unique(n, "rows");

len2 = 800

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
    plot3(ex_x,ex_y,ex_z,'o','markeredgecolor',ex_col,'markersize', 15, 'markerfacecolor', flip(c(ii,:)./255.0));
    ne_x = [p(ii,1)];
    ne_y = [p(ii,2)];
    ne_z = [p(ii,3)];
    plot3(ne_x,ne_y,ne_z,'o','markeredgecolor',nex_col,'markersize', 15, 'markerfacecolor', flip(p(ii,:)./255.0));
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
    plot3(ex_x,ex_y,ex_z,'o','markeredgecolor',ex_col,'markersize', 15, 'markerfacecolor', flip(c(ii,:)./255.0));
    blue = [blue; ex_x];
    redgreen = [redgreen; ex_y, ex_z];
end

plot3([0 255], [0 0], [0 0], 'b-');
plot3([0 0], [0 255], [0 0], 'g-');
plot3([0 0], [0 0], [0 255], 'r-');

xlabel ('Blue channel')
ylabel ('Green channel')
zlabel ('Red channel')

% Linear regression. Regress red/green against blue
b = blue\redgreen;
redgreen_calc = b .* blue;
plot3 (blue, redgreen_calc(:,1), redgreen_calc(:,2), 'k-')

% Add the fit line to fig 1 too
figure(1)
plot3 (blue, redgreen_calc(:,1), redgreen_calc(:,2), 'k-')