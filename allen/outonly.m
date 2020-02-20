b = unique(a, "rows");

len = size(b)(1)

len2 = 100
rndIDX = randperm(len);
c = b(rndIDX(1:len2), :);

figure (1)
clf
hold on
for ii = 1:len2
    ex_x = [0 c(ii,1)];%r
    ex_y = [0 c(ii,2)];%g
    ex_z = [0 c(ii,3)];%b
    plot3(ex_x,ex_y,ex_z,'ok-')
    %ne_x = [0 ne(ii,1)];
    %ne_y = [0 ne(ii,2)];
    %ne_z = [0 ne(ii,3)];
    %plot3(ne_x,ne_y,ne_z,'om-')
end

plot3([0 255], [0 0], [0 0], 'r-');
plot3([0 0], [0 255], [0 0], 'g-');
plot3([0 0], [0 0], [0 255], 'b-');
