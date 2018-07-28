function [mhat_lasso, lam_lasso, mhat_scad,lam_scad]=HDVCfit_qr(X,Y,U,N,p,hoptf, Cn, m,z,J,thre,iterbdd, tau)

%choosing the regularization parameter lamda_1 for the ordinary group lasso
%estimator by CV
folds =5;
splitsize = floor(N/folds)*ones(folds,1); %floor is a round function toward negative infinity
splitsize(folds) = N - (folds-1)*floor(N/folds); % I think it is identical with the above

lam_array=(0.075:0.025:0.125); % a seq 0.1 to 0.2 by 0.01
CVval_array=zeros(1,length(lam_array)); %
for lamk=1:length(lam_array) % grid search for lambda
rand('state', 1+m+lamk); % Initialize the seed for permutation sample
ptb = randperm(N); %random permutation for integer -> randomized sample
lam=lam_array(lamk);
temp=0;
txt = ['CV lamda is ', num2str(lam)];
disp(txt)

 for run=1:folds  % CV iteration
     txt = ['CV fold ', num2str(run), ' of ', num2str(folds)];
     disp(txt)
     tuneindex = ptb(((run - 1) * floor(N/folds) + 1) : ((run - 1) * floor(N/folds) + splitsize(run))); % the 1st cv index containing 40.
     idx_tu = zeros(N,1);
     idx_tu(tuneindex) = 1; % all zero but tune index is 1
        
     Ytu = Y(logical(idx_tu)); % logical is a function that assign non-zero value to logical 1 
     Utu = U(logical(idx_tu)); % and pick Z of which logical is 1.
     Xtu = X(logical(idx_tu),:); % separate data just for CV
        
     ztu=floor(Utu*J); % J is a discritization number 
     ztu(ztu==0)=1; % multiply 100 and when ztu==0, ztu=1 (Complementary)

     Ytr = Y(logical(1-idx_tu)); % not CV data
     Ntr = length(Ytr); 
     Xtr = X(logical(1-idx_tu),:);
     Utr = U(logical(1-idx_tu));
        
[mhattr, iter, error] = QPLVCM_LASSO(p, z, J,  Ntr, Xtr, Ytr, Utr, hoptf,lam, thre, iterbdd, tau);  %group lasso estimator        

mYtu=(Xtu*mhattr(1:p,:))+(Xtu*mhattr((1:p)+p,:)).*(outerop(Utu,-z,'+')/hoptf);   % N times J matrix (\hat Yi,-j (Zi:z))

mYt=zeros(length(Ytu),1);
for cvind=1:length(Ytu)
    mYt(cvind)=mYtu(cvind, ztu(cvind));
end

temp=temp+sum((Ytu-mYt).^2)/N;
 end
CVval_array(lamk)=temp;

end

txt = ['Final estimation in LASSO'];
disp(txt)

[optcv optind]=min(CVval_array); 

lam_lasso=lam_array(optind);
[mhat_lasso, iter, error] = QPLVCM_LASSO(p, z, J,  N, X, Y, U, hoptf,lam_lasso, thre, iterbdd, tau);  %group lasso estimator     

%choosing the regularization parameter lamda_1 for the ordinary group lasso
%estimator by BIC
lam_array=(0.075:0.025:0.125);
BIC_array=zeros(1,length(lam_array));
V_array=zeros(1,length(lam_array));
C_array=zeros(1,length(lam_array));

for lamk=1:length(lam_array) % select the best lambda with BIC
lam=lam_array(lamk);

txt = ['BIC lamda is ', num2str(lam)];
disp(txt)

% nu=[];
% w=[];
% for j=1:p
%     theta=norm(mhat0(j,:))/sqrt(J);
%     nu=[nu, Pprime(theta,lam)];
%     mconst=[mhat0(j,:)-mean(mhat0(j,:)), mhat0(j+p,:)];
%     %norm(mconst)/sqrt(J)
%     w=[w, Pprime(norm(mconst)/sqrt(J), lam)];
% end

[mhat, iter, error] = QPLVCM_SCAD(p, z, J,  N, X, Y, U, hoptf,lam, thre, iterbdd, tau); %jointly selection 


mY=(X*mhat(1:p,:))+(X*mhat((1:p)+p,:)).*(outerop(U,-z,'+')/hoptf);   
Yres=(Y*ones(1,J)-mY);
kw=norm_ker(outerop(U,z,'-'),hoptf);

V=nnz(mhat((1:p)+p,:).^2*ones(J,1));
C=nnz(mhat(1:p,:).^2*ones(J,1))-V;
V_array(lamk)=V;
C_array(lamk)=C;

BIC_array(lamk)=log(sum(sum(qr_obj((Yres.*kw), tau)))/(N*J))+(V*log(N*hoptf)/(N*hoptf)+C*log(N)/N)*Cn;

end


[optbic optind]=min(BIC_array); 


%final estimate
txt = ['Final estimation in SCAD'];
disp(txt)
lam_scad=lam_array(optind);

hlam2(m)=lam_scad; 

[mhat_scad, iter, error] = QPLVCM_SCAD(p, z, J,  N, X, Y, U, hoptf, lam_scad, thre, iterbdd, tau); %jointly selection 

end


