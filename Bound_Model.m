% model of distance
mu_X = 2.5;  % mu_G = 9
a = 2;     % \rho = 0.6
sigma_X = 0.6;  % sigma_G^2 = 1.8
gamma_th = 1;
f_egc = @(x_1,x_2) exp(a*x_1+x_2)+exp(x_1+a*x_2)-sqrt(2*gamma_th);
figure(1);
h = ezplot(f_egc,[-5,5,-5,5]);
C = get(h,'contourMatrix');
x_1 = C(1,2:end);
x_2 = C(2,2:end);

% upper bound
k = (a.*exp(a.*x_1+x_2)+exp(x_1+a.*x_2))./(exp(a*x_1+x_2)+a.*exp(x_1+a.*x_2));
D_upper =((mu_X-x_2)+k.*(mu_X-x_1))./sqrt(sigma_X^2.*(1+k.^2));
upper_bound = qfunc(D_upper);
figure;
subplot(1,2,1);
plot3(x_1,x_2,upper_bound);title('Upper bound');
xlabel('x_{m1}');
ylabel('x_{m2}');

% lower bound
lower_bound = qfunc((mu_X-x_1)./sigma_X).*qfunc((mu_X-x_2)./sigma_X);
subplot(1,2,2);
plot3(x_1,x_2,lower_bound);title('Lower bound');
xlabel('x_{m1}');
ylabel('x_{m2}');
