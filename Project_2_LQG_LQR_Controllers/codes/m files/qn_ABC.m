syms x(t) t1(t) t2(t) x_tt(t) t1_tt(t) t2_tt(t)
syms M m1 m2 l1 l2 g h F

x_t = diff(x,t);
 
t1_t = diff(t1,t);
 
t2_t = diff(t2,t);

% Using the Lagrangian computed in the report

L = (x_t^2*(M+m1+m2)/2) + (m1*(t1_t^2)*(l1^2))/2 + (m2*(t2_t^2)*(l2^2)/2) ...
    + (m1*l1*cos(t1)*(g-x_t*t1_t)) + (m2*l2*cos(t2)*(g-x_t*t2_t)) - M*g*h ;

eqn1 = diff(L,x) - diff(diff(L,x_t),t) + F == 0;

eqn2 = diff(L,t1) - diff(diff(L,t1_t),t) == 0;

eqn3 = diff(L,t2) - diff(diff(L,t2_t),t) == 0;

sol_t1_tt = isolate(eqn2, diff(t1,t,t));
sol_t2_tt = isolate(eqn3, diff(t2,t,t));

eqn1 = subs(eqn1, {diff(t1,t,t), diff(t2,t,t)}, {rhs(sol_t1_tt), rhs(sol_t2_tt)});

sol_x_tt = isolate(eqn1, diff(x,t,t));

sol_t1_tt = subs(sol_t1_tt, diff(x,t,t), rhs(sol_x_tt));

sol_t2_tt = subs(sol_t2_tt, diff(x,t,t), rhs(sol_x_tt));

A = sym('A%d%d', [6 6]);
B = sym('B%d%d', [6,1]);

A(1,1) = 0;
A(1,2) = 1;
A(1,3) = 0;
A(1,4) = 0;
A(1,5) = 0;
A(1,6) = 0;

A(2,1) = subs( diff(rhs(sol_x_tt), x), {x, t1, t2}, {0, 0, 0} );
A(2,2) = subs( diff(rhs(sol_x_tt), x_t), {x, t1, t2}, {0, 0, 0} );
A(2,3) = subs( diff(rhs(sol_x_tt), t1), {x, t1, t2}, {0, 0, 0} );
A(2,4) = subs( diff(rhs(sol_x_tt), t1_t), {x, t1, t2}, {0, 0, 0} );
A(2,5) = subs( diff(rhs(sol_x_tt), t2), {x, t1, t2}, {0, 0, 0} );
A(2,6) = subs( diff(rhs(sol_x_tt), t2_t), {x, t1, t2}, {0, 0, 0} );

A(3,1) = 0;
A(3,2) = 0;
A(3,3) = 0;
A(3,4) = 1;
A(3,5) = 0;
A(3,6) = 0;

A(4,1) = subs( diff(rhs(sol_t1_tt), x), {x, t1, t2}, {0, 0, 0} );
A(4,2) = subs( diff(rhs(sol_t1_tt), x_t), {x, t1, t2}, {0, 0, 0} );
A(4,3) = subs( diff(rhs(sol_t1_tt), t1), {x, t1, t2}, {0, 0, 0} );
A(4,4) = subs( diff(rhs(sol_t1_tt), t1_t), {x, t1, t2}, {0, 0, 0} );
A(4,5) = subs( diff(rhs(sol_t1_tt), t2), {x, t1, t2}, {0, 0, 0} );
A(4,6) = subs( diff(rhs(sol_t1_tt), t2_t), {x, t1, t2}, {0, 0, 0} );

A(5,1) = 0;
A(5,2) = 0;
A(5,3) = 0;
A(5,4) = 0;
A(5,5) = 0;
A(5,6) = 1;

A(6,1) = subs( diff(rhs(sol_t2_tt), x), {x, t1, t2}, {0, 0, 0} );
A(6,2) = subs( diff(rhs(sol_t2_tt), x_t), {x, t1, t2}, {0, 0, 0} );
A(6,3) = subs( diff(rhs(sol_t2_tt), t1), {x, t1, t2}, {0, 0, 0} );
A(6,4) = subs( diff(rhs(sol_t2_tt), t1_t), {x, t1, t2}, {0, 0, 0} );
A(6,5) = subs( diff(rhs(sol_t2_tt), t2), {x, t1, t2}, {0, 0, 0} );
A(6,6) = subs( diff(rhs(sol_t2_tt), t2_t), {x, t1, t2}, {0, 0, 0} );

B(1,1) = 0;
B(2,1) = subs( diff(rhs(sol_x_tt), F), {x, t1, t2}, {0, 0, 0} );
B(3,1) = 0;
B(4,1) = subs( diff(rhs(sol_t1_tt), F), {x, t1, t2}, {0, 0, 0} );
B(5,1) = 0;
B(6,1) = subs( diff(rhs(sol_t2_tt), F), {x, t1, t2}, {0, 0, 0} );

cont_mat = [B A*B A^2*B A^3*B A^4*B A^5*B];
determinant_cont = det(cont_mat);
