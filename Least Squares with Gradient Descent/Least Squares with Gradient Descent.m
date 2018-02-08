% 2. Least Squares with Gradient Descent:
clear all;
A = rand(3,2); %[7 3;5 8;6 2];
b = rand(3,1); %[3;1;4];
N_iter = 1;  % iteration index
tol = 1e-6;   % stopping criterion
x(:,1) = [0;0];  % starting point(from jy's work.how to compute starting point?)
t = 0.5;  % initial t

f = @(x) x'*A'*A*x-2*x'*A'*b+b'*b;   % objective function
grad_f = @(x) 2*A'*A*x-2*A'*b;       % gradient of the objective function

while( norm( grad_f(x(:,N_iter)) ) > tol )  
    % step 1: compute the direction "delta_x" at the current point:
            delta_x(:,N_iter) = -grad_f(x(:,N_iter));
    % step 2: find t through exact line search:
            t = (b-A*x(:,N_iter))./(-2)*A*A'*(A*x(:,N_iter)-b);
    %step 3: update x:
            x(:,N_iter+1) = x(:,N_iter) + t*delta_x(:,N_iter);
            N_iter = N_iter+1;
end
   
N_iter = N_iter-1;
% -------------- plotting contour ---------------
xx1 = min(x(1,1:N_iter)):0.01:max(x(1,1:N_iter));
xx2 = min(x(2,1:N_iter)):0.01:max(x(2,1:N_iter));
len1 = length(xx1);
len2 = length(xx2);
f_plot = zeros(len1,len2);
for i = 1:length(xx1);
    for j = 1:length(xx2)
        xx = [xx1(i);xx2(j)];
        f_plot(i,j) = f(xx);
    end
end
contour(f_plot,5);
xlabel('x(1)'),ylabel('x(2)');
grid on
hold on

% ------------ plotting iterating points ----------
for m = 1:N_iter
    plot(x(1,N_iter),x(2,N_iter),'o-');
    hold on
end

%     t = ( 2*delta_x(:,N_iter)'*A'*b - delta_x(:,N_iter)'*A'*A*x(:,N_iter) - x(:,N_iter)'*A'*A*delta_x(:,N_iter) )...
%             /2*delta_x(:,N_iter)'*A'*A*delta_x(:,N_iter);
