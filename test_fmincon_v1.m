clc
clear
close all hidden

% x0 = [1 0 0 1 1 0 1 1 1];
n = 9;
A = zeros(6,n);
A(1,[1 4]) = [1 -1];
A(2,[4 7]) = [1 -1];
A(3,[2 5]) = [1 -1];
A(4,[5 8]) = [1 -1];
A(5,[3 6]) = [1 -1];
A(6,[6 9]) = [1 -1];
b = zeros(6,1);

Aeq = zeros(3,n);
Aeq(1,1:3) = 1; Aeq(2,4:6) = 1; Aeq(3,7:9) = 1;
beq = [1;2;3];

lb = zeros(n,1)-0.1;         % lower bounds
ub = 1.1*ones(n,1);        % uppder bounds
fun = @(x) 1;   % objective function
nonlcon = @mycon;       % function handle to constraints
x0 = [1 0 0 1 1 0 1 1 1];             % initial guess
x0 = [1 1 1 1 1 1 1 1 1];             % initial guess

options = optimoptions('fmincon','Display','iter');     % display iterations
% call the solver
[x,fval,exitflag,~,lambda] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)


function [c,ceq] = mycon(x)
ceq(1) = sum(x(1:3).^2) - 1;
ceq(2) = sum(x(4:6).^2) - 2;
ceq(3) = sum(x(7:9).^2) - 3;
c = [];
end

