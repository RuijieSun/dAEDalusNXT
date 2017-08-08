%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [A0,A1,A2,E,D,gamma_new]=minimum_state_approximation_opt(Q,k,gamma,nits) 
%% rogers approximation
offset=0;

nQ=length(Q(:,1,1));
nL=length(gamma);
A_base=zeros(2*length(k)*nQ,3+nL);

Wr=ones(nQ,nQ,nL);
Wi=ones(nQ,nQ,nL);

%data normalization
epsi=1;
for i=1:nQ
    for j=1:nQ
        max_r = max(max(abs(real(Q(i,j,:)))),epsi);
        max_i = max(max(abs(imag(Q(i,j,:)))),epsi);
        Wr(i,j,:)=epsi/max_r;
        Wi(i,j,:)=epsi/max_i;
    end
end
%Weights added to remove the unrealistic coupling of the modes in approximation. Giving a
%high weight for the small values in Q forces the decoupling of the modes to be more accurately described. 
for i=1:nQ
    for j=1:nQ
        for k_idx=1:nL
            if abs(real(Q(i,j,k_idx)))<1E-4
                Wr(i,j,k_idx)=10;                                   
            end
            if abs(real(Q(i,j,k_idx)))<1E-5
                Wr(i,j,k_idx)=100;
            end
            if abs(real(Q(i,j,k_idx)))<1E-6
                Wr(i,j,k_idx)=1000;
            end
            if abs(real(Q(i,j,k_idx)))<1E-7
                Wr(i,j,k_idx)=10000;
            end
            if abs(real(Q(i,j,k_idx)))<1E-8
                Wr(i,j,k_idx)=100000;
            end 

            if abs(imag(Q(i,j,k_idx)))<1E-4
                Wi(i,j,k_idx)=10;
            end
            if abs(imag(Q(i,j,k_idx)))<1E-5
                Wi(i,j,k_idx)=100;
            end
            if abs(imag(Q(i,j,k_idx)))<1E-6
                Wi(i,j,k_idx)=1000;
            end
            if abs(imag(Q(i,j,k_idx)))<1E-7
                Wi(i,j,k_idx)=10000;
            end
            if abs(imag(Q(i,j,k_idx)))<1E-8
                Wi(i,j,k_idx)=100000;
            end 

        end
    end
end



E=zeros(nL,nQ);
D=ones(nQ,nL);
D(1:min(nL,nQ),1:min(nL,nQ))=eye(min(nL,nQ),min(nL,nQ));

    % Nonlinear Optimization (Sequential Simplex)
    options=optimset('Display','notify','MaxFunEvals',50,'TolX',1e-4);
    [gamma_new,fval,exitflag]=fminsearch(@cost,gamma,options);% fminsearch always minimize in respect to the first output
    exitflag
    gamma_new
    
    cost(gamma_new);
    
    function J=cost(gamma)

        for interations=1:1:nits
            redi=0;
            resnormi=0;
            for col_loop=1:1:nQ
                offset=0;
                for ii=1:1:nQ
                    k_idx=1;
                    for i=1:2:2*length(k)    
                        b(i+offset)=real(Q(ii,col_loop,k_idx))*Wr(ii,col_loop,k_idx);
                        b(i+1+offset)=imag(Q(ii,col_loop,k_idx))*Wi(ii,col_loop,k_idx);
                        A_base(i+offset,1+3*(ii-1):3*ii)=[1 0 -k(k_idx)^2]*Wr(ii,col_loop,k_idx);
                        A_base(i+1+offset,1+3*(ii-1):3*ii)=[0  k(k_idx) 0]*Wi(ii,col_loop,k_idx);
                        for j=1:1:nL
                            A_base(i+offset,3*nQ+j)=D(ii,j)*k(k_idx)^2/(k(k_idx)^2+gamma(j)^2)*Wr(ii,col_loop,k_idx);
                            A_base(i+offset+1,3*nQ+j)=D(ii,j)*k(k_idx)*gamma(j)/(k(k_idx)^2+gamma(j)^2)*Wi(ii,col_loop,k_idx);
                        end
                        k_idx=k_idx+1;
                    end
                    offset=offset+length(k)*2;
                end
                %% Solution for Rogers Approximation#
                [sol,nn,resid]=lsqlin(A_base,b');
                        redi=redi+sum(abs(resid));
                A0(:,col_loop)=sol(1:3:3*nQ);
                A1(:,col_loop)=sol(2:3:3*nQ);
                A2(:,col_loop)=sol(3:3:3*nQ);
                E(1:nL,col_loop)=sol(3*nQ+1:end);
            end
             redi   
            % D computation
            A_base=A_base*0;
            b=b*0;
            for row_loop=1:1:nQ
                offset=0;
                for ii=1:1:nQ
                    k_idx=1;
                    for i=1:2:2*length(k)
                        b(i+offset)=real(Q(row_loop,ii,k_idx))*Wr(row_loop,ii,k_idx);
                        b(i+1+offset)=imag(Q(row_loop,ii,k_idx))*Wi(row_loop,ii,k_idx);
                        A_base(i+offset,1+3*(ii-1):3*ii)=[1 0 -k(k_idx)^2]*Wr(row_loop,ii,k_idx);
                        A_base(i+offset+1,1+3*(ii-1):3*ii)=[0    k(k_idx) 0]*Wi(row_loop,ii,k_idx);
                        for j=1:1:nL
                            A_base(i+offset,3*nQ+j)=E(j,ii)*k(k_idx)^2/(k(k_idx)^2+gamma(j)^2)*Wr(row_loop,ii,k_idx);
                            A_base(i+offset+1,3*nQ+j)=E(j,ii)*k(k_idx)*gamma(j)/(k(k_idx)^2+gamma(j)^2)*Wi(row_loop,ii,k_idx);
                        end
                        k_idx=k_idx+1;
                    end
                    offset=offset+length(k)*2;
                end
                %% Solution for Rogers Approximation
                [sol,resnorm,resid]=lsqlin(A_base,b');
                    resnormi=resnormi+sum(resid.^2);
%                       resnormi=resnormi+sum((resid).^2./Wf(:,row_loop));
                A0(row_loop,:)=sol(1:3:3*nQ);
                A1(row_loop,:)=sol(2:3:3*nQ);
                A2(row_loop,:)=sol(3:3:3*nQ);
                D(row_loop,1:nL)=sol(3*nQ+1:end);
            end
            A_base=A_base*0;
            b=b*0;
            
        end
        
        J=sqrt(resnormi);
        
    end
    
    nQ=size(Q,1);
    debug=1;
    figure
    hold on
    if debug==1
        for i=1:6
            for j=7:16
                kk=1;
                for kc=min(k):0.01:max(k)
                    Qtilde(kk,1)=A0(i,j)-kc^2*A2(i,j);
                    Qtilde(kk,2)=kc*A1(i,j);
                    for lag=1:nL
                        Qtilde(kk,1)=Qtilde(kk,1)+E(lag,j)*D(i,lag)*kc^2/(kc^2+gamma_new(lag)^2);
                    end
                    
                    for lag=1:nL
                        Qtilde(kk,2)=Qtilde(kk,2)+E(lag,j)*D(i,lag)*kc*gamma_new(lag)/(kc^2+gamma_new(lag)^2);
                    end
                    kk=kk+1;
                end
                
                plot(real(squeeze(Q(i,j,:))),imag(squeeze(Q(i,j,:))),'-rx');
                
                plot(Qtilde(:,1),Qtilde(:,2),'-b');
                legend('Data','Rogers Approximation (N=6)')
            end
        end
    end
end