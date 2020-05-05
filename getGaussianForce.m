function f = getGaussianForce(A,s,x)
    f = (A./s.^2) .*x.*exp(-0.5.*(x./s).^2);