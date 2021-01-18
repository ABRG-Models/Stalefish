display('Test')
A=[ 0.1, 2.3, 1.8;
    0.8, 2.2, 1  ;
    1,   1,   1  ];

D=[ 0,1.9,2.1;
    1,2.1,1.2;
    1 , 1, 1];

# Let D = M*A

clf;
# Plot A red (the original data)
plot (A(1,:),A(2,:), 'ro-');
hold on;
# Plot D blue - the transformed data
plot (D(1,:), D(2,:), 'bo-')

# Now some example vectors
# 1.89, 0.98 to 2.08, 0.91
p1=[1.89,0.98,1];
p2=[2.08,0.91,1];
plot ([1.89, 2.08], [0.98, 0.91], 'r-')
# 2.1 1.55 to 2.32, 1.5
p3=[2.1,1.55,1];
p4=[2.32,1.5,1];
plot ([2.1, 2.32], [1.55, 1.5], 'r-')
# 2.34,2.17 to 2.54,2.173
p5=[2.34,2.17,1];
p6=[2.54,2.173,1];
plot ([2.34, 2.54], [2.17, 2.173], 'r-')

p7=[1.107,1.49,1]
p8=[0.85,1.67,1]
plot ([p7(1), p8(1)], [p7(2), p8(2)], 'r-')

# How to transform point p2 to p2': p2p = p2 * (A\D); p2p = p2p ./ p2p(3);
M=D * inv(A)

p1p = M * p1';
p2p = M * p2';
plot ([p1p(1),p2p(1)],[p1p(2),p2p(2)],'b')

p3p = M * p3';
p4p = M * p4';
plot ([p3p(1),p4p(1)],[p3p(2),p4p(2)],'b')

p5p = M * p5';
p6p = M * p6';
plot ([p5p(1),p6p(1)],[p5p(2),p6p(2)],'b')

p7p = M * p7';
p8p = M * p8';
plot ([p7p(1),p8p(1)],[p7p(2),p8p(2)],'b')
