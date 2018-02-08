% project 1:compressive sensing
clear all;
n = 150;  % 50,100,500  n:length of vector x
S = 3;  %  S<<n  S:number of non-zeros in x 
k = 5;  %  k:length of vector y
error = zeros(1,n);

% generating x:
x0 = zeros(n,1);  
S_index = ceil( rand(S,1)*n );  % determine the indexes of non-zeros 
x0(S_index) = rand()*n;

for k = 1:n
    % generating phi:
        phi = rand(k,n);
        for i = 1:k
            for j = 1:n
                if phi(i,j) > 0.5
                    phi(i,j) = 1;
                else
                    phi(i,j) = -1;
                end
            end
        end

        y = phi * x0;      % take the measurement

        cvx_begin
            variable x(n);
            minimize norm(x,1);
            subject to
            phi*x == y;
        cvx_end
        
        error(1,k) = sum( abs(x-x0) );
end

figure;
plot(error);
grid on
title('Error between recovered x and original x (n=150)');
xlabel('k');
ylabel('error');



% results:
% figure(1);
% plot(abs(x-x0),'r');
% title('Error between recovered x and original x');
% legend('error','original x0','recovered x');
% hold on
% grid on
% plot(x0,'g');
% plot(x,'b');
