%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
%
% The code originates mainly from Proteus
% (https://github.com/sbind/JointTool/tree/master/ext/prOOteus) uploaded by
% Mario Natella

classdef class_bucklingElement
    %CLASS_BUCKLINGELEMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %dimension in x-coordinate of buckling element coordinate sys
        elLength
        %dimension in y-coordinate of buckling element coordinate sys
        elWidth
        %constant K terms
        stiff
        %number of critical modes
        nModes=2
        %specialType (e.g. shearOnly)
        type='';
        %lamId 
        lamId=0;
    end
    
    methods
        %constructor
        function obj=class_bucklingElement(l,w)
            % calc constant K terms Matrix
            obj.elLength=l;
            obj.elWidth=w;
            obj=obj.calcConstantStiff();   
        end
        
        function obj=calcConstantStiff(obj)
             x=[0;obj.elLength;obj.elLength;0];
             y=[0;0;obj.elWidth;obj.elWidth];
             n=10;%Lobatto Parameter
             [L,L1,L2,wg,sg] = obj.shapefun(n);
             ng = length(sg);

             K11 = zeros(n*(n+1)/2); K12 = K11; K16 = K11; K22 = K11; K26 = K11; K66 = K11;
             Kxx = K11; Kxy = K11; Kyy = K11;

             for iy = 1:ng
                 for ix = 1:ng
                     [xg,Ag,Jig,Hgx,Hgy] = obj.mapping(x,y,sg(ix),sg(iy));
                     w = zeros(n*(n+1)/2,1);
                     thx = w; thy = w; kxx = w; kxy = w; kyy = w; 
                     k = 0;
                     for r = 1:n
                         for t = 1:r
                            k=k+1;
                            Pval = L(t,ix)*L(r-t+1,iy);
                            Psx = L1(t,ix)*L(r-t+1,iy); Psy = L(t,ix)*L1(r-t+1,iy);
                            Psxx = L2(t,ix)*L(r-t+1,iy); Psxy = L1(t,ix)*L1(r-t+1,iy); Psyy = L(t,ix)*L2(r-t+1,iy);
                            [w(k),thx(k),thy(k),kxx(k),kxy(k),kyy(k)] = obj.tansformvalues(Pval,Psx,Psy,Psxx,Psxy,Psyy,Jig,Hgx,Hgy);
                         end
                     end
                     f = wg(ix)*wg(iy)*Ag;
                     K11 = K11 + kxx*kxx'*f;
                     K12 = K12 + (kxx*kyy'+kyy*kxx')*f;
                     K16 = K16 + (kxx*kxy'+kxy*kxx')*f;
                     K22 = K22 + kyy*kyy'*f;
                     K26 = K26 + (kyy*kxy'+kxy*kyy')*f;
                     K66 = K66 + kxy*kxy'*f;
                     Kxx = Kxx - thx*thx'*f;
                     Kxy = Kxy - (thx*thy'+thy*thx')*f;
                     Kyy = Kyy - thy*thy'*f;
                 end
             end

             obj.stiff.K11 = K11;
             obj.stiff.K12 = K12;
             obj.stiff.K16 = K16;
             obj.stiff.K22 = K22;
             obj.stiff.K26 = K26;
             obj.stiff.K66 = K66;
             obj.stiff.Kxx = Kxx;
             obj.stiff.Kxy = Kxy;
             obj.stiff.Kyy = Kyy;
            
            
            
        end
        
        
        function r=calcBucklingReserveFactors(obj,strains,ABD)
            if strcmp(obj.type, 'shearOnly')
                strains(1:2,:)=strains(1:2,:)*0;
            end
            N=ABD(1:3,1:3)*strains;
            %loop over all strains
            
            r=zeros(obj.nModes,size(strains,2));
            for iStrain=1:size(strains,2)
                [r(:,iStrain)]=obj.quadbuckl(N(:,iStrain),ABD(4:6,4:6));             
            end
        end
        function [r] = quadbuckl(obj,N,D)
            nmodes=obj.nModes;
            K11 = obj.stiff.K11;
            K12 = obj.stiff.K12;
            K16 = obj.stiff.K16;
            K22 = obj.stiff.K22;
            K26 = obj.stiff.K26;
            K66 = obj.stiff.K66;
            Kxx = obj.stiff.Kxx;
            Kxy = obj.stiff.Kxy;
            Kyy = obj.stiff.Kyy;

            % Phi = zeros(3,3,nmodes);
            Kg = N(1)*Kxx+N(2)*Kyy+N(3)*Kxy;
            K = D(1,1)*K11+D(1,2)*K12+D(1,3)*K16+D(2,2)*K22+D(2,3)*K26+D(3,3)*K66;


            [~,r] = obj.poseig(K,Kg,nmodes);
        end  
        
    end
    
    methods(Static)
        function [L,L1,L2,w,x] = shapefun(n)

            if n < 1
                L = []; L1 = []; L2 = []; w = []; x = [];
                return
            end


            % Gauss weights and abscissas
            T = diag([n:-1:1]./sqrt(2*[n+1:-1:2]-1)./sqrt(2*[n:-1:1]-1),-1); T = T + T';
            [V,xg] = eig(T);
            w = V(end,:).^2;
            x = diag(xg)';
            w = 2*w/sum(w);
            % Lobatto shape functions and first and second derivatives
            L = zeros(n,length(x)); L1 = L; L2 = L;

            L2(1,:) = 1;
            L1(1,:) = x;
            L(1,:) = (x.*L1(1,:)-1)/2;

            if n == 1, return; end

            L2(2,:) = 3*x;
            L1(2,:) = 3*L(1,:)+1;
            L(2,:) = (x.*L1(2,:)-L1(1,:))/3;

            for k = 3:n
                L1(k,:) = ((2*k-1)*x.*L1(k-1,:)-(k-1)*L1(k-2,:))/k;
                L(k,:) = (x.*L1(k,:)-L1(k-1,:))/(k+1);
                L2(k,:) = k*L1(k-1,:)+x.*L2(k-1,:);
            end

        end
        
        function [xg,Ag,Jig,Hgx,Hgy] = mapping(x,y,sx,sy)

            Ng = [(1-sx)*(1-sy)/4,(1+sx)*(1-sy)/4,(1+sx)*(1+sy)/4,(1-sx)*(1+sy)/4];
            N1g = [-(1-sy)/4,(1-sy)/4,(1+sy)/4,-(1+sy)/4;
                -(1-sx)/4,-(1+sx)/4,(1+sx)/4,(1-sx)/4];            
            xg = Ng*x; yg = Ng*y;
            Jg = [N1g*x N1g*y];
            Hgx = ([1/4,-1/4,1/4,-1/4]*x)*[0 1;1 0];
            Hgy = ([1/4,-1/4,1/4,-1/4]*y)*[0 1;1 0];         
            Ag = Jg(1,1)*Jg(2,2)-Jg(2,1)*Jg(1,2);
            Jig = inv(Jg);

        end
        
        function [w,thx,thy,kxx,kxy,kyy] = tansformvalues(Pval,Psx,Psy,Psxx,Psxy,Psyy,Jig,Hx,Hy)

            w = Pval;
            thx = Jig(1,1)*Psx+Jig(1,2)*Psy;
            thy = Jig(2,1)*Psx+Jig(2,2)*Psy;
            H = Jig*([Psxx Psxy;Psxy Psyy]-thx*Hx-thy*Hy)*Jig';
            kxx = H(1,1);
            kyy = H(2,2);
            kxy = H(1,2)+H(2,1);

        end 
        
        function [pVk,pek] = poseig(K,Kg,k)

        n = size(K,1);
            [V,D] = eig(Kg,K);
            e = diag(D); % eigenvalues
            pei = find(e>0); % indices of positive eigenvalues
            pn = length(pei); % number of positive eigenvalues
            pe = e(pei);
            [pe,pes] = sort(-pe); % sorting
            pei = pei(pes);

            if(k<=pn)
                pek = -pe(1:k);
                pVk = V(:,pei(1:k));
            else
                pek(1:pn) = -pe(1:pn);
                pek(pn+1:k) = 0;
                pVk = V(:,pei(1:pn));
                pVk(1:n,pn+1:k) = 0;
            end
        end 
    end
end

