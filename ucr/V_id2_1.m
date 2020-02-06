load V_Id2_1.h5

s = [Frame000.sbox_centers Frame001.sbox_centers Frame002.sbox_centers Frame003.sbox_centers ...
     Frame004.sbox_centers Frame005.sbox_centers Frame006.sbox_centers Frame007.sbox_centers ...
     Frame008.sbox_centers Frame009.sbox_centers Frame010.sbox_centers Frame011.sbox_centers ...
     Frame012.sbox_centers Frame013.sbox_centers Frame014.sbox_centers Frame015.sbox_centers ...
     Frame016.sbox_centers Frame017.sbox_centers Frame018.sbox_centers ]';

m = [Frame000.means Frame001.means Frame002.means Frame003.means ...
     Frame004.means Frame005.means Frame006.means Frame007.means ...
     Frame008.means Frame009.means Frame010.means Frame011.means ...
     Frame012.means Frame013.means Frame014.means Frame015.means ...
     Frame016.means Frame017.means Frame018.means ]';

figure (1)
scatter3 (s(:,1),s(:,2),s(:,3),400*m,m,"filled")

figure (2)
clf;
hold on;
xx = Frame000.class.layer_x * ones (1,size(Frame000.sbox_linear_distance)(2));
scatter (xx, Frame000.sbox_linear_distance,  100, Frame000.means, "filled")

xx = Frame001.class.layer_x * ones (1,size(Frame001.sbox_linear_distance)(2));
scatter (xx, Frame001.sbox_linear_distance,  100, Frame001.means, "filled")

xx = Frame002.class.layer_x * ones (1,size(Frame002.sbox_linear_distance)(2));
scatter (xx, Frame002.sbox_linear_distance,  100, Frame002.means, "filled")

xx = Frame003.class.layer_x * ones (1,size(Frame003.sbox_linear_distance)(2));
scatter (xx, Frame003.sbox_linear_distance,  100, Frame003.means, "filled")

xx = Frame004.class.layer_x * ones (1,size(Frame004.sbox_linear_distance)(2));
scatter (xx, Frame004.sbox_linear_distance,  100, Frame004.means, "filled")

xx = Frame005.class.layer_x * ones (1,size(Frame005.sbox_linear_distance)(2));
scatter (xx, Frame005.sbox_linear_distance,  100, Frame005.means, "filled")

xx = Frame006.class.layer_x * ones (1,size(Frame006.sbox_linear_distance)(2));
scatter (xx, Frame006.sbox_linear_distance,  100, Frame006.means, "filled")

xx = Frame007.class.layer_x * ones (1,size(Frame007.sbox_linear_distance)(2));
scatter (xx, Frame007.sbox_linear_distance,  100, Frame007.means, "filled")
