% Separating Two Sets of Points via Feasibility of an LP
% this piece of code is under the help of Yangjun LI (My teammate). Thanks a lot!

% Initialization
X = [-1,-0.8,-3;1,-0.5,-2];
Y = [5,4;9,-2];
a = [1,1]';
b = 1;
alpha = 0.1;
beta = 0.5;
mu = 2;
e = 10^(-3);
t = 1;

s = select_largest(X,Y,a,b);
g = gradient_sep(s,a,b,X,Y,t);
x = [s,a',b]';
iteration1 = 0;
iteration2 = 0;

while s>0
    iteration1 = iteration1 + 1;
    iteration2 = 0;
    while g'*g>e
        t1 = 1;
        iteration2 = iteration2 + 1;
        iter2(iteration1) = iteration2;
        while f_sep(s-t1*g(1),a-t1*g(2:3),b-t1*g(4),X,Y,t)>f_sep(s,a,b,X,Y,t)-alpha*t1*g'*g
            t1 = beta*t1;
        end
        x = x - t1*g;
        x_value(:,iteration2)=x;
        g = gradient_sep(s,a,b,X,Y,t);
        g_value(:,iteration2)=g;
    end
    s = x(1);
    a(1) = x(2);
    a(2) = x(3);
    b = x(4);
    t = mu*t;
    g = gradient_sep(s,a,b,X,Y,t);
    optimal_s(iteration) = s;
end
