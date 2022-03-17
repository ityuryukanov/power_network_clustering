function [MK, K, M] = svm_em(pst)
% 
% Compute the linearized electromechanical model.
% All loads are modeled as constant impedance. 
% 
jay = sqrt(-1);

% Save & set up global variables
pststate = PST.savepst();
clearvars -global;
PST.pst_var;
bus = pst.bus(:,1:10);   %remove extra-added data to avoid bugs
line = pst.lin(:,1:min(size(pst.lin,2),10));   %remove extra-added data to avoid bugs
mac_con = pst.gen;
if size(mac_con,2)<20
  assert(norm(mac_con(:,[4:6,8:15,18]),'fro')<1e-6,'Classical generator models are required');  
else
  assert(norm(mac_con(:,[4:6,8:15,18,20:21]),'fro')<1e-6,'Classical generator models are required');    
end
basmva = pst.basmva;
basrad = 2*pi*pst.f0; 
%syn_ref = 1;   %machine 1 (arbitrary) is reference
[num_mach, dummy] = size(mac_con);

%Solve input loadflow
tol  = 1e-9;    % tolerance for convergence
iter_max = 300; % maximum number of iterations
vmin = 0.5;     % voltage minimum
vmax = 1.5;     % voltage maximum
acc  = 1.0;     % acceleration factor
[DUMMY,bus,lin_flw] = evalc('PST.LLoadflow(bus,line,tol,iter_max,vmin,vmax,acc,''n'',2);');

% step 1: construct reduced Y matrix - all loads are assumed
%         to be of constant impedance load
Y_red = PST.RED_YYbus(bus,line);     % pre-fault
f = PST.mac_em(0,1,bus,0);          % computes globals for mac_ang, eqprime, edprime..!

assert(all(abs(edprime)<1e-4));
eprime = abs(edprime+jay*eqprime); % E' of generators
Kb = zeros(num_mach,num_mach);     %synchronizning torque coefficients Part 1
Kg = zeros(num_mach,num_mach);     %synchronizning torque coefficients Part 2
for i=1:1:num_mach
  for j=1:1:num_mach
    if i~=j
      Bij = imag(Y_red(i,j));
      Gij = real(Y_red(i,j));
      ang = mac_ang(i) - mac_ang(j);
      Kb(i,j) = eprime(i)*eprime(j)*Bij*cos(ang);
      Kg(i,j) = eprime(i)*eprime(j)*(-Gij*sin(ang));
    end
  end
  Kb(i,i) = -sum(Kb(i,:));
  Kg(i,i) = -sum(Kg(i,:));  
end
assert(isreal(Kg));
assert(isreal(Kb));
%%Kb = (Kb+Kb')/2;                   %ensure symmetry??
K = Kb + Kg;
H_rated = mac_con(:,16);
S_rated = mac_con(:,3);
H = H_rated.*S_rated/basmva;       %to the common power base
M = sparse(diag(2*H/basrad));      %M=2Hi/w0
MK = M\K;                          %MK=inv(M)*K;
assert(isreal(MK));
K = cat(3, K, Kb, Kg);

%Restore the original PST state
PST.loadpst(pststate);
end

