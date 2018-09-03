data = dlmread('Q_dist_output.dat', '', 1, 0);
numNaN = sum(sum(isnan(data)))/3;

Ql = 1.09999999973479;
Q0 = 1;
alpha = 0.01;

(Ql-Q0)^2/alpha
