

          % e_z   % e_h   % e_mi  % e_me  % e_a  % e_l    % e_qk   % e_y% %% e_bh     % e_be   %e_d    %e_r_t    %e_eps_K_b 
Std_e = [0.027     0      0       0       0       0        0        0      0         0        0          0          0              % e_z
          0      0.076    0       0       0       0        0        0      0         0        0          0          0              % e_h
          0        0     0.003    0       0       0        0        0      0         0        0          0          0              % e_mi
          0        0      0      0.007    0       0        0        0      0         0        0          0          0              % e_me
          0        0      0       0      0.06     0        0        0      0         0        0          0          0              % e_a
          0        0      0       0       0     0.57       0        0      0         0        0          0          0              % e_l
          0        0      0       0       0       0       0.019     0      0         0        0          0          0              % e_qk
          0        0      0       0       0       0        0      0.634    0         0        0          0          0              % e_y
          0        0      0       0       0       0        0        0     0.067      0        0          0          0              % e_bh 
          0        0      0       0       0       0        0        0      0         0.063    0          0          0              % e_be
          0        0      0       0       0       0        0        0      0         0        0.033      0          0              %e_d
          0        0      0       0       0       0        0        0      0         0        0          0.01      0              %e_r_t 0.002
          0        0      0       0       0       0        0        0      0         0        0          0          0.031];        %e_eps_K_b
                

shocks;
var e_ae  =  Std_e(5,5)^(2);
var wbh   =  0.05^2;
var e_z   =  Std_e(1,1)^(2);
var e_h   =  Std_e(2,2)^(2);
var e_mi  =  Std_e(3,3)^(2);
var e_me  =  Std_e(4,4)^(2);
var e_l   =  Std_e(6,6)^(2);
var e_qk  =  Std_e(7,7)^(2);
var e_y   =  Std_e(8,8)^(2);
var e_bh  =  Std_e(9,9)^(2);
var e_be  =  Std_e(10,10)^(2);
var e_d   =  Std_e(11,11)^(2);
var e_r_t =  Std_e(12,12)^(2);
var e_eps_K_b =Std_e(13,13)^(2);
end;

stoch_simul(order=1, irf = 40, nograph);
//r_t r_be r_bh pi l_i l_p y c inv d_p spr_b K_b;