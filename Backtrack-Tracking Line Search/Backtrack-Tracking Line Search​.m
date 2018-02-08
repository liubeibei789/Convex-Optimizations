% plot objective function
clear all;
xx1 = -0.7:0.02:0.3;
xx2 = -0.3:0.02:0.3;
[X1,X2]=meshgrid(xx1,xx2);
f_plot = exp(X1+2*X2) + exp(X1-2*X2) + exp(-X1);
mesh(X1,X2,f_plot);
contour(X1,X2,f_plot,5);
xlabel('x1'),ylabel('x2');
title('exp(X1+2*X2) + exp(X1-2*X2) + exp(-X1)');
grid on
hold on

f = @(x) exp(x(1)+2*x(2)) + exp(x(1)-2*x(2)) + exp(-x(1));         % objective function
grad_f = @(x) [exp(x(1)+2*x(2)) + exp(x(1)-2*x(2)) - exp(-x(1)) ;  % gradient of the objective function
               2*exp(x(1)+2*x(2)) - 2*exp(x(1)-2*x(2)) ];
hessian_f = @(x) [exp(x(1)+2*x(2)) + exp(x(1)-2*x(2)) + exp(-x(1)) 2*exp(x(1)+2*x(2)) - 2*exp(x(1)-2*x(2)) ;
                  2*exp(x(1)+2*x(2)) - 2*exp(x(1)-2*x(2)) 4*exp(x(1)+2*x(2)) + 4*exp(x(1)-2*x(2)) ];
                                                                   % hessian of the objective function
alpha = 0.1;
beta = 0.5;
x(:,1) = [0.5;0.5];
t = 1;
N_iter = 1;
tol = 1e-6;    % termination of algorithm

while( norm( 0.5*grad_f(x(:,N_iter))'*power(hessian_f(x(:,N_iter)),-1)*grad_f(x(:,N_iter)) ) > tol )
    % step 1:compute delta_x:
            delta_x(:,N_iter) = -power(hessian_f(x(:,N_iter)),-1)*grad_f(x(:,N_iter));
    % step 2:compute t through backtracking line search:
            t = 1;
            while( f( x(:,N_iter)+t*delta_x(:,N_iter) ) >= f(x(:,N_iter)) + alpha*t*grad_f(x(:,N_iter))'*delta_x(:,N_iter) )
                t = t*beta;
            end
    % step 3: update x:
            x(:,N_iter+1) = x(:,N_iter) + t*delta_x(:,N_iter);
            N_iter = N_iter+1;
end

N_iter = N_iter-1;
% ------------ plotting iterating points ----------
for m = 1:N_iter
    plot(x,'o-');
    hold on
end
