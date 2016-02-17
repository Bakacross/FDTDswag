% Constant (void)
u = 4*pi*1e-12;
eps = 8.8541878176e-12;

% Initial condition

hy = [1,1];
ez = [0,1];

for i = 1:100
    hy(i+0.5,i) = hy1d(i,i);
end