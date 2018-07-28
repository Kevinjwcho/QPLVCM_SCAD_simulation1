function [mhat, iter, error] = QPLVCM_LASSO(p, z, J,  N, X, Y, U, hoptf,lam, thre, iterbdd, tau )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% containing varying ceofficient values for a and b
mhat=zeros(2*p,J);   % (m_1(z)^T;...;m_p(z)^T;m_1^(1)(z)^T;...;m_p^(1)(z)^T)

iter=1;
while(1==1) % I think it is for break; iterating everlast
    
omhat=mhat;    % zero matrix is omhat
%%%%%%%%%%%%%%%%%%%%%%
for j=1:p % estimate parameters columnwise

    %%%%Generate Yres first
    mY=(X*mhat(1:p,:))+(X*mhat((1:p)+p,:)).*(outerop(U,-z,'+')/hoptf);   % Erase hoptf
    mYj=X(:,j)*mhat(j,:)+(X(:,j)*mhat(j+p,:)).*(outerop(U,-z,'+')/hoptf);  
    Yres=(Y*ones(1,J)-(mY-mYj)).*norm_ker(outerop(z,U,'-'), hoptf)';  % N times J matrix (Yi-\hat Yi,-j (Zi:z)) I think mYj is not necessary but I don't know the solution because of 'for' sentence
    Yres=reshape(Yres', 1, N*J); % rowise reshape
    
    %%% Xj part %%%%
    Xj_full = zeros(2*J, J*N);
    for grd_N=1:N
        Xj_row = zeros(2*J, J); % initial
        for grd_J=1:J
            slope = ((U-z(grd_J))/hoptf);
            ker = norm_ker(z(grd_J)-U, hoptf);
            Xv=[X(grd_N,j) X(grd_N,j)*slope(grd_N)]*ker(grd_N);
            in=Xv';
            Xj_row((1:2)+(grd_J-1)*2, grd_J)= in;   
        end
        Xj_full(:, (1:J)+(grd_N-1)*J) = Xj_row;
    end

T = [eye(2*J),zeros(2*J,1)]; %eye: identity matrix
T1 = ones(N*J,1); %eye: identity matrix
f = [zeros(1,2*J),lam];

A = [eye(2*J)/sqrt(J), zeros(2*J,1)]; % basis vectors
V= [zeros(1,2*J), 1];


cvx_begin
     variable B(2*J+1);
     minimize ((1/(N*J))*(qr_obj((Yres-(B'*T')*Xj_full), tau)*T1)+f*B)
       subject to 
       norm(A*B) <= V*B;
  %     eta_p >= 0;
   %    eta_n >= 0;
cvx_end

if (B(2*J+1)<=10^-6) % to make simple
    B=zeros(2*J+1,1);
end

for grd_J=1:J
mhat(j,grd_J)=B(2*(grd_J)-1); %checking the respective norm of (m_j, m_j^{(1)}), j=1,...,p.
mhat(j+p,grd_J)=B(2*(grd_J));
end

end
error=sum(sum((omhat-mhat).^2))/J;   % stop if error is \sum_{j=1}^p \int (m_j^{(k)}-m_j^{(k-1)})^2 <=p *10^{-5}

   if(error<thre)
      disp('error is')                 
      disp(error)                   
    break;
   end
   
   txt = [num2str(iter), 'th iters in QPLVCM_LASSO in processing'];
disp(txt)
iter=iter+1;

if(iter>=iterbdd)                 
      disp('error is')                 
      disp(error)                     
      break;
end

end



end

