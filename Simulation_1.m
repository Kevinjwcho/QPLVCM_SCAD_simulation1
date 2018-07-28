cvx_setup
cvx_quiet true
%Choosing the bandwidth based on the test sample with N=5000
tau=0.3;
 p=7;
 p0=p-1;
 Co = 0.5.^abs(outerop(1:p0,1:p0,'-')); %covariance matrix
 
thre=10^-5*p;
iterbdd=20;

 J=100; %discritization number
 z=(1/J:1/J:1); %sequence from 0 to 1 by 0.01


hoptf=0.15; % bandwidth problem

% Data settings
 N = 300;
 M= 100; % # of whole analysis iteration

% ind_nz0=zeros(M,p);
% ISE0=zeros(M,p);
% ind_nz=zeros(M,p);
% ind_v=zeros(M,p);
% ISE=zeros(M,p);
% 
% hlam1=zeros(M,1); 
% hlam2=zeros(M,1);

ind_nz0_qr=zeros(M,p);
ind_v0_qr=zeros(M,p);
ISE0_qr=zeros(M,p);
ind_nz_qr=zeros(M,p);
ind_v_qr=zeros(M,p);
ISE_qr=zeros(M,p);

hlam1_qr=zeros(M,1); 
hlam2_qr=zeros(M,1);


% Cn=1;
Cn=1/2;
%Cn=sqrt(log(p));
%Cn=log(p);

tic %start stopwatch timer
for m=1:M
txt = ['Whole iter ', num2str(m), ' of ', num2str(M)];
disp(txt);
rand('state', 1+m); % Initialize the seed for the covariates (fix seed)
randn('state', 1+m); %Initialize the seed for the covariates (fix seed)

data = Data_GN_S1(N,Co,p0);

U = data{2};
X = data{3};
Y=data{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mhat0_qr, lam0_qr, mhat_qr,lam_qr]=HDVCfit_qr(X,Y,U,N,p,hoptf,Cn,m,z,J,thre,iterbdd,tau);
 % mhat0 represents varying coeffs with only lam 1
ind_nz0_qr(m,:)=1-(mhat0_qr(1:p,J)==0); 
ind_v0_qr(m,:)=1-(mhat0_qr(p+1:p+p,J)==0);   % # of significant varying coefficients
ISE0_qr(m,1)=sum((mhat0_qr(1,:)-a1(z)).^2)/J;  % MSE
ISE0_qr(m,2)=sum((mhat0_qr(2,:)-a2(z)).^2)/J;
ISE0_qr(m,3)=sum((mhat0_qr(3,:)-a3(z)).^2)/J;
ISE0_qr(m,4)=sum((mhat0_qr(4,:)-a4(z)).^2)/J;
ISE0_qr(m,5)=sum((mhat0_qr(5,:)-a5(z)).^2)/J;

for k=6:p
   ISE0_qr(m,k)=(sum(mhat0_qr(k,:).^2))/J; 
end

hlam1_qr(m)=lam0_qr;
ISE_qr(m,1)=sum((mhat_qr(1,:)-a1(z)).^2)/J;
ISE_qr(m,2)=sum((mhat_qr(2,:)-a2(z)).^2)/J;
ISE_qr(m,3)=sum((mhat_qr(3,:)-a3(z)).^2)/J;
ISE_qr(m,4)=sum((mhat_qr(4,:)-a4(z)).^2)/J;
ISE_qr(m,5)=sum((mhat_qr(5,:)-a5(z)).^2)/J;

for k=6:p
   ISE_qr(m,k)=(sum(mhat_qr(k,:).^2))/J;
end

ind_nz_qr(m,:)=1-(mhat_qr(1:p,J)==0);      % # of significant nonzero coefficients
ind_v_qr(m,:)=1-(mhat_qr(p+1:p+p,J)==0);   % # of significant varying coefficients

hlam2_qr(m)=lam_qr;


end
toc
save Simulation_1.mat
%Saving the results

csvwrite('result/ISE0_qr.csv',ISE0_qr);
csvwrite('result/ind_nz0_qr.csv',ind_nz0_qr);
csvwrite('result/ind_v0_qr.csv',ind_v0_qr);
csvwrite('result/hlam0_qr.csv',hlam1_qr);

csvwrite('result/ISE_qr.csv',ISE_qr);
csvwrite('result/ind_nz_qr.csv',ind_nz_qr);
csvwrite('result/ind_v_qr.csv',ind_v_qr);
csvwrite('result/hlam2_qr.csv',hlam2_qr);




