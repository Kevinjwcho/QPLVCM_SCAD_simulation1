function y=Pprime(theta,lambda)

a=3.7;
y=(theta<=lambda)*lambda + (a*lambda - theta)/(a-1) .*((theta>lambda)&(theta<a*lambda)) + 0*(theta>a*lambda);



