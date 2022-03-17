% file: pst_var.m
%
% Syntax: pst_var
%
% Purpose: Define global variables for power system simulation. 
%         

%
% Input:
%
% Output:
%
% Files:
%
% See Also:

% Algorithm:
%
% Calls:
%
% Call By:

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved

% History (in reverse chronological order)
%
% Version:
% Date:
% Author:
% Purpose:
% Modification:

% Version:  1.0
% Author:   Joe H. Chow
% Date:     January 1991

% system variables
global  basmva basrad syn_ref mach_ref sys_freq
global  bus_v bus_ang psi_re psi_im cur_re cur_im bus_int

% synchronous machine variables
global  mac_con mac_pot mac_int
global  mac_ang mac_spd eqprime edprime psikd psikq
global  curd curq curdg curqg fldcur
global  psidpp psiqpp vex eterm theta ed eq 
global  pmech pelect qelect
global  dmac_ang dmac_spd deqprime dedprime dpsikd dpsikq

% excitation system variables
global  exc_con exc_pot Efd V_R V_A V_As R_f V_FB V_TR V_B
global  dEfd dV_R dV_As dR_f dV_TR
global  exc_sig

% non-conforming load variables
global  load_con load_pot load_int

% svc variables
global  svc_con svc_pot B_cv dB_cv
global  svc_sig

% pss variables
global  pss_con pss_pot pssvar3
global  pss1 pss2 pss3 dpss1 dpss2 dpss3

% turbine-governor variables
global  tg_con tg_pot 
global  tg1 tg2 tg3 dtg1 dtg2 dtg3 

% tcsc variables
global  tcsc_con tcsc_pot tcsc_sig
global  tcsc1 tcsc2 tcsc3 dtcsc1 dtcsc2 dtcsc3
global  tcsc_in

% lin1 variables
global lin1_con lin1_pot var4
global lin1a lin1b lin1c lin1d dlin1a dlin1b dlin1c dlin1d

% lin2 variables
global lin2_con lin2_pot lin2c_ou
global lin2_1a lin2_2a lin2_1b lin2_2b lin2_1c lin2_2c
global dlin2_1a dlin2_2a dlin2_1b dlin2_2b dlin2_1c dlin2_2c

% DC variables?
global dcr_dc dci_dc svc_dc tcsc_dc;
