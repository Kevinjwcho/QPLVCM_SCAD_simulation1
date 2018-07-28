function f=norm_ker(u,h)
%f=normpdf(x,0,h);
f=(1-(u/h).^2).*(abs(u/h)<=1)*3/(4*h);

end