% generate 100 iid obserations from normal distribution, set mu as 7
n = 1000;
error = random('norm',0,1,[1,n-1]);
X0 = zeros(1,n);
X0(1) = 2;
a = 1; b = 0.6;

for i = 2:n
    X0(i) = a + b*X0(i-1) +error(i-1);
end
Xs = [X0(1:n-1);X0(2:n)];

Lfunction = @(X,theta) ((1/sqrt(2*pi))* exp(-((X(1,:)-(theta(2)+theta(1)*X(2,:))).^2)/2));

[X,FVAL] = MLE(Xs,[1,1],Lfunction) 

J = [1 mean(X0(1:n-1)); mean(X0(1:n-1)) mean(X0(1:n-1).^2)];
Lambda = inv(J)

Z_alpha = 1.96;  
Interval = [X' - sqrt(diag(Lambda))/sqrt(n)*Z_alpha, X' + sqrt(diag(Lambda))/sqrt(n)*Z_alpha]

% P-value: here we use Chi-Squared test
% H0: theta_hat = theta_0
theta_0 = [0.6 1];
Eta = n*(X-theta_0)*inv(Lambda)*(X-theta_0)';
p = 1-chi2cdf(Eta,2)
