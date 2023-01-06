syms m1 g m2 M L1 L2
m1 = 100;
m2 = 100;
M = 1000;
L1 = 20;
L2 = 10;
g = 9.81;
q0 = [2 0 deg2rad(30) 0 deg2rad(10) 0];
tspan = 0:0.1:100;

A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -((M*g)+(m1*g))/(M*L1) 0 -g*m2/(M*L1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*L2) 0 -((M*g)+(m2*g))/(M*L2) 0];
B = [0; 1/M; 0; 1/(L1*M); 0; 1/(L2*M)];
C1 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];
C2 = [0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];
C3 = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 1 0];
C4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];

D = [0; 0; 0];

Obs1_rank = rank([C1' A'*C1' ((A')^2)*C1' ((A')^3)*C1' ((A')^4)*C1' ((A')^5)*C1']);
Obs2_rank = rank([C2' A'*C2' ((A')^2)*C2' ((A')^3)*C2' ((A')^4)*C2' ((A')^5)*C2']);
Obs3_rank = rank([C3' A'*C3' ((A')^2)*C3' ((A')^3)*C3' ((A')^4)*C3' ((A')^5)*C3']);
Obs4_rank = rank([C4' A'*C4' ((A')^2)*C4' ((A')^3)*C4' ((A')^4)*C4' ((A')^5)*C4']);

sys1 = ss(A,B,C1,D);
sys3 = ss(A,B,C3,D);
sys4 = ss(A,B,C4,D);

fprintf('Rank of observability matrix for Output vector 1 (x) = %d' , Obs1_rank);

fprintf('Rank of observability matrix for Output vector 2 (theta1, theta2) = %d' , Obs2_rank);

fprintf('Rank of observability matrix for Output vector 3 (x, theta2) = %d' , Obs3_rank);

fprintf('Rank of observability matrix for Output vector 4 (x, theta1, theta2) = %d' , Obs4_rank);