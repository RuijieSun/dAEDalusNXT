%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Authors: Klaus Seywald (klaus.seywald@mytum.de) 
%          and Simon Binder (simon.binder@tum.de)
% 
% This file is part of dAEDalusNXT (https://github.com/seyk86/dAEDalusNXT)
% Partly taken from https://github.com/sbind/JointTool/tree/master/ext/prOOteus/structures 
%
function obj = cross_sectional_modeler(obj)
% constants
E = [0 0 -1 0;0 1 0 0];
Ie = zeros(4,6);
Ie(1,1)= 1;Ie(2,5)= 1;Ie(3,6)= 1;Ie(4,4)= 1;
Is = zeros(2,6);
Is(1,2)= 1;Is(2,3)= 1;


% Nr of dof per node (only 4 relevant unknown shell warping displacements,
% see finite element discretization appendix of Abdalla's cross sectional
% modeler paper for more info)
ndof = 4;

% Nr of shell elements within the cross section
Ne = size(obj.elements,1);

% Nr of nodes within the cross section
Nn = size(obj.nodes_yz,1);

% calculates the total nr of DOFs within the cross section, i.e. total
% amount of unique nodes times ndof
dof_max = ndof*Nn;




% runs assemble function. This function initializes the matrices and
% data necessary to run the actual modeller
[F,H0,K00,H1,K10,B0,G] = assemble(Ne,dof_max,obj);
% Looping over all elements, see comments above for Element idx meaning



% creates arrays with constrained DOFs and free DOFs
dof_tot       = 1:dof_max; %total dof

% Clamping a whole node in order to avoid singular stiffness matrix.
node_clamped  = 1; %constrained node ID 
dof_clamped   = (node_clamped-1)*ndof+1:ndof*node_clamped;%constrained dof array
dof_free      = setdiff(dof_tot,dof_clamped); % free dof array


% first order approximation of warping w0 and Euler stiffness matrix
V0 = zeros(size(H0));
V0(dof_free,:) = -K00(dof_free,dof_free)\H0(dof_free,:);
%euler bernoulli beam elements
obj.Se = F+H0'*V0; % Euler: modulus/stiffness matrix
obj.Se = (obj.Se+obj.Se')/2; %Symmetrize
% second order approximation of warping w1 and the Timoshenko stiffness matrix
% (for the moment these are not used, but they will turn out to be useful if we upgrade the
% Daedalus beam model from an Euler Bernoulli one to a Timoshenko one).
V1 = zeros(size(H0));
H1b = H1+K10*V0;
V1(dof_free,:) = K00(dof_free,dof_free)\(H1b(dof_free,:));

P=H0'*V1;
obj.Ce =  obj.Se\eye(4);                     %Euler: compliance/flexibility matrix
obj.Cs =  E*obj.Ce*(V1'*H1b+P'*obj.Ce*P)*obj.Ce*E';  %Shear: compliance of the shear stiffness components
obj.Ces = obj.Ce*P*obj.Ce*E';                    %Coupling: euler-shear forces

C6 = Ie'*obj.Ce*Ie+Is'*obj.Cs*Is-(Ie'*obj.Ces*Is+Is'*obj.Ces'*Ie); % Timoshenko: compliance/flexibility matrix

S6 = C6\eye(6);  % Timoshenko: modulus/stiffness matrix
obj.S6 = (S6+S6')/2; %Symmetrize
%% Normalized strain variation over cross-section: Gamma = Gamma_hat*epsilon

% This is matrix used to recover the shell strains from the Nodal deflections
obj.Gamma_euler = (G+B0*V0);

% This is matrix used to recover the shell strains from the reaction forces
Gamma_euler2 = (G+B0*V0)*obj.Ce;

% The 2 matrices below are not used for now. If we upgrade to a Timoshenko model, they will be needed
obj.Gamma_shear = (B0*V1-Gamma_euler2*P)*obj.Ce;
obj.Gamma       = Gamma_euler2*Ie+obj.Gamma_shear*E'*Is; %strain across cross-seaction

end

% This function initializes the matrices and data necessary to run the actual modeller
function [F,H0,K00,H1,K10,B0,G] = assemble(Ne,dof,obj)


K00 = zeros(dof,dof);
H0  = zeros(dof,4);
H1  = zeros(dof,4);
K10 = zeros(dof,dof);
F   = zeros(4,4);
R   = zeros(dof,6); % maxtrix where the columns are the 6 ridig body modes of a cross-section

B0  = zeros(6*Ne,dof); %B0 matrix of all elements
G   = zeros(6*Ne,4); %G matrix for all elements
L   = 0;

% looping through all shell elements within the cross section
for el = 1:Ne
    node1Id=obj.elements(el,2);
    node2Id=obj.elements(el,3);
    elmdof = [(node1Id-1)*4+1:node1Id*4 (node2Id-1)*4+1:node2Id*4 ];  % DOFs associated with nodes of shell element "el"
    nelm = 6*(el-1)+1:6*el;  % Indices of strains associated with current element (DOFs of element itselfm useful for strain recovery)
 

    % yz coords (cross sectional modeller coord sys) of both nodes of
    % current shell element
    yz_el_loc = obj.nodes_yz([node1Id, node2Id],:);
    
    % if statements below determine if current element belongs to fs,
    % sk_up, rs or sk_lo. Depending on the result, the corresponding
    % ABD_stiff stiffness matrix is assigned to the el_ABD variable. The
    % shell_id_arr corresponding to the current element is also saved within an array.
    if obj.elements(el,4)==1
        el_ABD = obj.laminate_fs.ABD_stiff;
    elseif obj.elements(el,4)==2
        el_ABD = obj.laminate_sk_up.ABD_stiff;
    elseif obj.elements(el,4)==3
        el_ABD = obj.laminate_rs.ABD_stiff;
    elseif obj.elements(el,4)==4
        el_ABD = obj.laminate_sk_lo.ABD_stiff;
    elseif obj.elements(el,4)==5
        el_ABD = obj.laminate_st.ABD_stiff;
    end
    
    %run the function that calculates local properties, assembling all
    %output matrices of the cross sectional modeller.
    [Fe,H0e,K00e,H1e,K10e,Re,B0e,Ge,Le] = elmvar_reduced(yz_el_loc,el_ABD);
    
    %assemble H,K,etc.....
    H0(elmdof,:)       = H0(elmdof,:)+H0e;
    K00(elmdof,elmdof) = K00(elmdof,elmdof)+K00e;
    H1(elmdof,:)       = H1(elmdof,:)+H1e;
    K10(elmdof,elmdof) = K10(elmdof,elmdof)+K10e;
    F                  = F+Fe;
    L                  = L+Le;
    %strain comp
    R(elmdof,:)        = R(elmdof,:)+Re;
    B0(nelm,elmdof)    = B0e;
    G(nelm,:)          = Ge;
end
end


function [F,H0,K00,H1,K10,R,B0,G,L] = elmvar_reduced(yz_coords,el_ABD_in)
% this function calulates the coefficients of element strain energy

%% element properties
el_ABD = el_ABD_in;          %material property of an element
y = yz_coords(:,1); %y coordinate of the 2 nodes
z = yz_coords(:,2); %z coordinate of the 2 nodes

L = sqrt(diff(y)^2+diff(z)^2); % Length of the shell element
ydot = diff(y)/L; %y derivative wrt s
zdot = diff(z)/L; %z derivative wrt s
%% the rigid body modes of currect element: ri i1,2,..6 for 8dof
r1 = [1,0,0,0,1,0,0,0]';
r2 = [0,ydot,-zdot,0,0,ydot,-zdot,0]';
r3 = [0,zdot,ydot,0,0,zdot,ydot,0]';
r4 = [0,y(1)*zdot-z(1)*ydot,z(1)*zdot+y(1)*ydot,-2,0,y(2)*zdot-z(2)*ydot,z(2)*zdot+y(2)*ydot,-2]';
r5 = [z(1),0,0,0,z(2),0,0,0]';
r6 = -[y(1),0,0,0,y(2),0,0,0]';

%% B0 and B1 and G
B00 = [0 0 0 0 0 0 0 0;
    0 -1/L 0 0 0 1/L 0 0;
    -1/L 0 0 0 1/L 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 -6/L^2 4/L 0 0 6/L^2 2/L;
    0 0 0 0 0 0 0 0];

B01 = [0 0 0 0 0 0 0 0;
    0 -1/L 0 0 0 1/L 0 0;
    -1/L 0 0 0 1/L 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 6/L^2 -2/L 0 0 -6/L^2 -4/L;
    0 0 0 0 0 0 0 0];

B10 = [1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 -3/(2*L) -1/4 0 3/(2*L) 3/4]; 

B11 = [0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 -3/(2*L) 3/4 0 0 3/(2*L) -1/4];


%% normal notation of the rigid body modes
G0 = [1 z(1)  -y(1)      0;
    0 0     0         0;
    0 0     0 (y(1)*zdot-z(1)*ydot);
    0 ydot  zdot     0;
    0 0     0         0;
    0 0     0        -2];
G1 = [1 z(2)  -y(2)      0;
    0 0     0         0;
    0 0     0 (y(2)*zdot-z(2)*ydot);
    0 ydot  zdot     0;
    0 0     0         0;
    0 0     0        -2];

%% transformation matrix from local to global

% if id_elem == 1
%     T4 = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
%     T = [T4 zeros(4); zeros(4) T4];
% elseif id_elem == 2
%     T4 = [1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1];
%     T = [T4 zeros(4); zeros(4) T4];
% elseif id_elem == 3
%     T4 = [1 0 0 0; 0 0 -1 0; 0 -1 0 0; 0 0 0 1];
%     T = [T4 zeros(4); zeros(4) T4];
% elseif id_elem == 4
%     T4 = [1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 1];
%     T = [T4 zeros(4); zeros(4) T4];
% end

T = eye(8);
T(2,2) = ydot;
T(2,3) = -zdot;
T(3,2) = zdot;
T(3,3) = ydot;
T(6,6) = ydot;
T(6,7) = -zdot;
T(7,6) = zdot;
T(7,7) = ydot;


%% coefficients of the strain energy expression
F   = L*(   1/3*G0'*el_ABD*G0    +    1/6*(G0'*el_ABD*G1  +  G1'*el_ABD*G0)+    1/3*G1'*el_ABD*G1);
H0  = L*(    1/3*B00'*el_ABD*G0     +     1/6*B00'*el_ABD*G1   +   1/6*B01'*el_ABD*G0   +   1/3*B01'*el_ABD*G1);
K00 = L*(     1/3*B00'*el_ABD*B00  +  1/6*(B00'*el_ABD*B01  +  B01'*el_ABD*B00)  +   1/3*B01'*el_ABD*B01);
H1  = L*(  1/3*B10'*el_ABD*G0  +  1/6*B10'*el_ABD*G1  +  1/6*B11'*el_ABD*G0  +  1/3*B11'*el_ABD*G1);
K10 = L*(1/3*B10'*el_ABD*B00+1/6*(B10'*el_ABD*B01+B11'*el_ABD*B00)+1/3*B11'*el_ABD*B01);
    %K11 =L*(1/3*B1i'*C*B1i+1/6*(B1i'*C*B1j+B1j'*C*B1i)+1/3*B1j'*C*B1j);
    
    %% transformation to the global(cartesian) coordinate system
    H0 = T*H0;
    K00 = T*K00*T';
    H1 = T*H1;
    K10 = T*K10*T';
    
    %% some additional parameters
    R = T*[r1 r2 r3 r4 r5 r6]; %matrix of rigid body modes(in global coordinate system)
    B0 = 1/2*(B00+B01)*T'; %geometry info strain contribution of warping displacement
    B1 = (1/2)*(B11+B10);
%     G = B1*[r1 r5 r6 r4]; %1/2*(G0+G1); %strain contribution of the beam strains
    G = (1/2)*(G1+G0);
end