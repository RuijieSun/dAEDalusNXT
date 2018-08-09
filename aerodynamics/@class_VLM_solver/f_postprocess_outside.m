%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
function [ F_body,M_body, F_aero, cp, cl , M_BA] = f_postprocess_outside( Uinf, Ma_corr, panels, Cind, Gamma, grid, rho, qinf, S_ref )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% as second input either 'quick' for no
            alpha=atan(Uinf(3)/Uinf(1));
            beta=atan(Uinf(2)/Uinf(1));
            if Ma_corr<1
                beta_inf=Ma_corr;
            else
                beta_inf=1;
            end
            M_BA=[  cos(alpha)*cos(beta) , -cos(alpha)*sin(beta), -sin(alpha);
                sin(beta)             , cos(beta)            ,    0;
                sin(alpha)*cos(beta) , -sin(alpha)*sin(beta), cos(alpha);];
            
      
			%determine downwash at the one quarter point
            nPanels = size(panels,2);
            u=(Cind(:,1:nPanels)*Gamma)/(4*pi);
            v=Cind(:,nPanels+1:2*nPanels)*Gamma/(4*pi);
            w=Cind(:,2*nPanels+1:3*nPanels)*Gamma/(4*pi);
            F_body=zeros(3,nPanels);
            F_aero=zeros(3,nPanels);
            M_body=zeros(3,nPanels);
            cp=zeros(1,nPanels);
            cl=zeros(3,nPanels);
            for i=1:1:length(Gamma)

                p_i=0.75*grid(:,panels(1,i))+0.25*grid(:,panels(4,i));
                p_o=0.75*grid(:,panels(2,i))+0.25*grid(:,panels(3,i));
                sp=p_o-p_i;

                %obj.crossp(:,i)=(rho*Gamma(i)*cross(sp',Uinf'))';
                
%                 Uinf(3)=Uinf(3)+uvwind(3); 
%                 Uinf(2)=Uinf(2)+uvwind(2);
%                 Uinf(1)=Uinf(1)+uvwind(1);
%                 
                % do this without using cross
                F_body(:,i)=1/beta_inf*(rho*Gamma(i)*cross(sp',(Uinf+[u(i),v(i),w(i)])'))';
                % like this:
%                 F_body(:,i)=1/beta_inf*(rho*Gamma(i)*[sp(2)*(Uinf(3)+w(i))-sp(3)*(Uinf(2)+v(i))  sp(3)*(Uinf(1)+u(i))-sp(1)*(Uinf(3)+w(i)) sp(1)*(Uinf(2)+v(i))-sp(2)*(Uinf(1)+u(i));])';
                %nastran way of force computation
%                 F_body(:,i)=1/beta_inf*(rho*Gamma(i)*norm(sp)*obj.colloc_nvec(:,i)'*norm((Uinf+0*[u(i),v(i),w(i)])))';

                M_body(:,i)=zeros(3,1);
               % F_body(2,i)=-F_body(2,i);
                
                %alpha2=atan((Uinf(3)+uvwind(3))/(Uinf(1)+uvwind(1)));
                %beta=atan((Uinf(2)+uvwind(2))/(Uinf(1)+uvwind(1)));
                
                F_aero(:,i)=M_BA'*F_body(:,i);
                % F_aero(:,i)=(rho*Gamma(i)*cross(sp',Uinf'))';
                % obj.cdi(:,i)=-(rho*Gamma(i)*obj.wind(i)*cross(sp',Uinf'/norm(Uinf)))';
                % obj.cdi(:,i)=M_BA'*(rho*Gamma(i)*cross(sp',(Uinf'+uvwind)/norm(Uinf+uvwind')))';
                sc=(norm(grid(:,panels(4,i))-grid(:,panels(1,i)))+norm(grid(:,panels(3,i))-grid(:,panels(2,i))))/2;
                cp(i)=1/beta_inf*2*Gamma(i)/(norm(Uinf)*sc);
                
                %obj.A(i)=norm(cross(grid(:,panels(2,i))-grid(:,panels(1,i)),grid(:,panels(4,i))-grid(:,panels(1,i))));
                cl(:,i)=F_body(:,i)/(qinf*S_ref);
            end
            
      

end

