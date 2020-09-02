function func1 = likelihood(alpha0, beta1,alpha1,Y)
      sigma0 = (1/length(Y)*sum(Y.^2));
      sigma = sqrt(alpha0 + beta1*Y(1)^2 + alpha1*sigma0.^2);
      func1 = log(1/sqrt(2*pi)*sigma)*exp(-Y(1).^2/(2*sigma.^2));
      t = length(Y);
      
      for i = t
          sigma = sqrt(alpha0 + beta1*Y(i-1).^2 + alpha1*sigma.^2);
          func1 = func1 + log(1/(sqrt(2*pi)*sigma)*exp(-Y(i).^2/(2*sigma.^2)));
      end
end