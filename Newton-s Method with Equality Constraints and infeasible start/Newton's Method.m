% Newton's Method with Equality Constraints and infeasible start

% Newton's Method with Equality Constraints 
% (Implement the infeasible start Newton’s method)
clear all;
A = [1 -1];
b = [0;0];

f = @(x) exp(x(1)+2*x(2)) + exp(x(1)-2*x(2)) + exp(-x(1));          % objective function
grad_f = @(x) [exp(x(1)+2*x(2)) + exp(x(1)-2*x(2)) - exp(-x(1)) ;   % gradient of the objective function
               2*exp(x(1)+2*x(2)) - 2*exp(x(1)-2*x(2)) ];
hessian_f = @(x) [exp(x(1)+2*x(2)) + exp(x(1)-2*x(2)) + exp(-x(1)) 2*exp(x(1)+2*x(2)) - 2*exp(x(1)-2*x(2)) ;
                  2*exp(x(1)+2*x(2)) - 2*exp(x(1)-2*x(2)) 4*exp(x(1)+2*x(2)) + 4*exp(x(1)-2*x(2)) ];   
                                                                    % hessian of the objective function
r_dual = @(x,mu) grad_f(x)+A'*mu;        % residual function: dual part
r_pri = @(x,mu) A*x-b;                   % residual function: primal part

x(:,1) = [0.5;0.5];  % starting point
alpha = 0.1;   
beta = 0.5;
mu = 1;
t = 1;
N_iter = 1;     % iteration index
tol = 1e-4;     % termination of algorithm

while( norm( [r_dual(x(:,N_iter),mu(:,N_iter))  r_pri(x(:,N_iter),mu(:,N_iter))] ) > tol || A*x(:,N_iter)~= b )
    % step 1: compute delta_x and delta_mu:
            delta_x(:,N_iter) = -power(hessian_f(x(:,N_iter)),-1) * ( grad_f(x(:,N_iter))+A'* delta_mu(:,N_iter) );
            delta_mu(:,N_iter) = power(( A* power(hessian_f(x(:,N_iter)),-1) *A'),-1) * ( -b+A*x(:,N_iter)-A*power(hessian_f(x(:,N_iter)),-1)*grad_f(x(:,N_iter)) );
    % step 2: find t through finding the norm of residual:   
            while(  norm( [ r_dual( x(:,N_iter)+t*delta_x(:,N_iter), mu(:,N_iter)+t*delta_mu(:,N_iter) )  r_pri(x(:,N_iter),mu(:,N_iter))] ) > tol )
                t = t*beta;
            end
    % step 3: update x and mu:
            x(:,N_iter+1) = x(:,N_iter) + t*delta_x(:,N_iter);
            mu(:,N_iter+1) = mu(:,N_iter) + t*delta_mu(:,N_iter);
            N_iter = N_iter+1;
end

 
% ------------ plotting iterating points ----------
for m = 1:N_iter
    plot(x(:,N_iter),'o-');
    hold on
end
