%%%QKLASSO settings(1)
function dat = Data_GN_S1(n,Co, p)


Mu = zeros(1,p);
%x = mvnrnd(Mu,Co,n);
x = [ones(n,1) mvnrnd(Mu,Co,n)];
% e = trnd(3,n,1);
% gm=gmdistribution([0;0], cat(3, 1, 25), [0.8, 0.2]);
% e = random(gm, n);
e = normrnd(0,1,n,1);
u = unifrnd(0,1,n,1);
y = a1(u).*x(:,1) + a2(u).*x(:,2)+a3(u).*x(:,3)+a4(u).*x(:,4)+a5(u).*x(:,5)+e;

dat{1} = y;
dat{2} = u;
dat{3} = x;

end