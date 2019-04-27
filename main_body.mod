

var kb_rat, wB, wbhk,
@#if amort
rho_amk,
@#endif
lambda_p $\lambda^{p}$, c_p $c^{p}$, h_p $h^{p}$, q_h $q^{p}$, d_p $d^{p}$, w_p $W^{p}$, l_p $l^{p}$                              %patient
lambda_i $\lambda^{i}$, c_i $c^{i}$, h_i $h^{i}$, s_i $s^{i}$, w_i $W^{i}$, l_i $l^{i}$, b_i $b^{i}$   
lambda_e $\lambda^{e}$, c_e $c^{e}$, q_k $q^{k}$, s_e $s^{e}$, y_e $y^{e}$, r_k $r^{k}$, u, k_e $k^{e}$, x, l_ep $l^{ep}$, l_ei $l^{ei}$, b_ee $b^{ee}$      %entreprenuers 
pi $\pi$, pi_wp $\pi^{wp}$, pi_wi $\pi^{wi}$, PIW                                               %
epsilon_d $\epsilon^{d}$, epsilon_be $\epsilon^{be}$, epsilon_bh $\epsilon^{bh}$, epsilon_qk $\epsilon^{qk}$, epsilon_l $\epsilon^{l}$, epsilon_y $\epsilon^{y}$, epsilon_z $\epsilon^{z}$, epsilon_h $\epsilon^{h}$, m_e $m^{e}$, a_e $a^{e}$, m_i $m^{i}$, eps_K_b $\varepsilon^{Kb}$  
k, inv, c, y, j_r $j^{r}$, y1, c_g $c^{g}$, tax, c_pr $c^{pr}$,                                % market clearing  
R_b $R^{b}$,r_bh $r^{bh}$, r_be $r^{be}$, r_t $r^{t}$, r_d $r^{d}$, K_b $Kb$, B, BH, BE, D, b_h $b^{h}$, b_e $b^{e}$, d_b $d^{b}$, j_b $j^{b}$, J_b $J^{b}$, v_b $v^{b}$, spr_b $spr^{b}$, weight_BH $w^{BH}$, weight_BE $w^{BE}$,  % banking
y_obs, c_obs, i_obs, q_h_obs, r_t_obs, r_be_obs, r_bh_obs, r_d_obs, b_e_obs, b_h_obs, d_b_obs, pi_obs, w_i_obs, nega;  
varexo vb_ss, m_i_ss, wbh, rho_am, e_z $e^{z}$, e_h $e^{h}$,  e_mi $e^{mi}$,  e_me $e^{me}$,  e_ae $e^{ae}$,  e_l $e^{l}$,   e_qk $e^{qk}$,  e_d $e^{d}$,  e_y $e^{y}$,  e_bh $e^{bh}$,  e_be $e^{be}$,  e_r_t $e^{rt}$,  e_eps_K_b $e^{epsKb}$;

parameters 
myltv
hab_p  $a^{p}$
beta_p  $\beta^{p}$
hab_i $a^{i}$
beta_i $\beta^{i}$  
hab_e  $a^{e}$
delta  $\delta$
csi1  $\chi^{1}$
csi2  $\chi^{2}$
myalp  $\myalp$
mu  $\mu$
beta_e   $\beta^{e}$
k_i   $k^{i}$
iota_w   $\iota^{w}$
piss  $\pi^{ss}$
k_be  $k^{be}$
phi $\phi$
k_p $k^{p}$
k_w $k^{w}$
k_kb $k^{kb}$
k_bh  $k^{kh}$
k_d $k^{d}$
gamma_p $\gamma^{p}$
gamma_i $\gamma^{i}$
gamma_e $\gamma^{e}$
gamma_b  $\gamma^{b}$
phi_R  $\phi^{R}$
phi_pi $\phi^{pi}$  
phi_y $\phi^{y}$
delta_b $\delta^{b}$
iota_p $\iota^{p}$
rho_e_z  $\rho^{ez}$
rho_e_h  $\rho^{eh}$
rho_mi  $\rho^{mi}$
rho_me  $\rho^{me}$
rho_e_ae $\rho^{eae}$
rho_e_l  $\rho^{el}$
rho_e_qk  $\rho^{eqk}$
rho_e_y  $\rho^{ey}$
rho_e_bh  $\rho^{ebh}$
rho_e_be  $\rho^{ebe}$
rho_e_d   $\rho^{ed}$
rho_eps_K_b $\rho^{epsKb}$
//m_i_ss   $m^{iss}$
m_e_ss  $m^{ess}$
epsilon_d_ss  $\epsilon^{dss}$
epsilon_be_ss $\epsilon^{bess}$ 
epsilon_bh_ss  $\epsilon^{bhss}$
epsilon_y_ss   $\epsilon^{yss}$
epsilon_l_ss  $\epsilon^{lss}$
h 
eps_d  $\epsilon^{d}$
eps_bh $\epsilon^{bh}$
eps_be $\epsilon^{be}$
pi_share $\pi^{share}$
r_t_ss  $r^{tss}$
r_d_ss  $r^{dss}$
r_bh_ss  $r^{bhss}$
r_be_ss  $r^{bess}$
r_k_ss   $r^{kss}$
J  
eps_b $\epsilon^{b}$
//vb_ss $vb^{ss}$ 
eta $\eta$
pcg2pcp;

myltv = 0;
hab_p  =0;
beta_p =0.99455; 
hab_i  =0.99;0
beta_i =0.975;  
hab_e  =0;
delta =0.03838; 
myalp =0.34; 
mu  = 0.8;
beta_e = beta_i;
phi = 1.5;
//vb_ss  = 0.09*1.02;
eta    = 0.20*1.250;
iota_w =0.2757; 
iota_p =0.1605 ;           
k_p   =28.6502;

@#if kbfree
  k_kb  = 15;0;
@#else
  k_kb  =11.0683;
@#endif

@#if kldfree
  k_bh  =0;
  k_be  =0;
  k_d   =0;
@#else
  k_bh  =10.0867;
  k_be  =9.3638;
  k_d   =3.5030;
@#endif

k_i   =10.1822;
k_w   =99.8983;
gamma_p =1;
gamma_i =1;
gamma_e =1;
gamma_b =1;
phi_R  =0.75;
phi_pi  =1.9816;
phi_y   =0.3459;           
rho_e_z  =0.95;
rho_e_h  =0.95;
rho_mi  =0.95;
rho_me  =0.95;
rho_e_ae =0.95;
rho_e_l  =0;
rho_e_qk  =0.95;
rho_e_y  =0;
rho_e_bh  =0;
rho_e_be  =0;
rho_e_d =0;
rho_eps_K_b  = 0.90;
piss  = 1;                                                         %steady state gross inflation rate
eps_d = - 1.1;                                                     % elast. of subst. of deposits 
eps_be = 5.85;                                                     % elast. of subst. of loans to H
eps_bh=  4.85;                                                     % elast. of subst. of loans to E

@#if freebank
  epsilon_d_ss     = 2.71828;                            % steady state  mark down interest rate on deposit   =0.5935
  epsilon_be_ss    = 2.71828;                           % steady state mark up on loans to E    =1.4717
  epsilon_bh_ss    = 2.71828;                           % steady state mark up on loans to H    =1.5587
@#else
  epsilon_d_ss     = (eps_d/(eps_d - 1));                            % steady state  mark down interest rate on deposit   =0.5935
  epsilon_be_ss    =(eps_be/(eps_be - 1));                           % steady state mark up on loans to E    =1.4717
  epsilon_bh_ss    =(eps_bh/(eps_bh - 1));                           % steady state mark up on loans to H    =1.5587
@#endif


//m_i_ss =0.75;                                                      % steady state loan-to-value ratio impatient households
m_e_ss =0.40;                                                      % steady state loan-to-value ratio Entrepreneurs

epsilon_y_ss = 8.9;
epsilon_l_ss = 5; 
r_d_ss  = (piss/beta_p) -1 ;
r_t_ss  = (piss/beta_p -1) *(eps_d - 1)/eps_d;                     % steady state gross nominal interest rate
r_bh_ss = epsilon_bh*r_t_ss ;                                      % steady state interest rate on loans to E
r_be_ss = epsilon_be*r_t_ss ;                                      % steady state interest rate on loans to H
r_k_ss  = -(1-delta)-0.4*(1-delta)*piss/beta_e*
               (1/(1+r_be_ss)-beta_e/piss)+1/beta_e;               % steady state rental rate of capital
eps_b        = 0.5*(eps_bh + eps_be);
pi_share    = 0.55;                   % share of bank profits payed out to patient households OK!!
delta_b = r_t_ss/0.09* (eps_d - eps_b + v_b*eps_d*(eps_b-1))/((eps_b-1)*(eps_d-1)); 
csi1  = r_k_ss;
csi2  =0.1*r_k_ss;
h =1;                                                              % fixed supply housing
J = 0.203;
pcg2pcp = 0.32; 

model;
nega =  e_ae;
exp(lambda_p) = (1 -hab_p) * exp(epsilon_z) / (exp(c_p) - hab_p * exp(c_p(-1)));                                %%%/* dL / d c_p */
J*(exp(epsilon_h)) / exp(h_p)  - exp(lambda_p) * exp(q_h) + beta_p * exp(lambda_p(+1)) * exp(q_h(+1)) =0;       %%%/* dL / d h_p */
exp(lambda_p) = beta_p *  exp(lambda_p(+1)) * (1 + exp(r_d))/exp(pi(+1));                                       %%%/* dL / d d_p */
exp(c_p) + exp(q_h) * ( exp(h_p) - exp(h_p(-1))) + exp(d_p) = exp(w_p) * exp(l_p) + (1 + exp(r_d(-1)))* exp(d_p(-1))/exp(pi) + exp(j_r) / gamma_p + (1-pi_share)*exp(j_b(-1))/(gamma_p*exp(pi));         %%%/* BC*/
exp(lambda_i) = (1 - hab_i) * exp(epsilon_z) / ( exp( c_i) - hab_i * exp( c_i(-1))) ;                           %%%/* dL / d c_i */
J*( exp(epsilon_h)) / exp(h_i) - exp(lambda_i) * exp(q_h) + beta_i * exp(lambda_i(+1)) * exp(q_h(+1)) + exp(s_i) * exp( m_i) * exp(q_h(+1)) * exp(pi(+1)) = 0 ;     %%%  /* dL / d h_i*/
exp(lambda_i)  -  beta_i *  exp( lambda_i ( +1)) *( 1 + exp ( r_bh)) / exp(pi(+1)) = exp(s_i) * ( 1 + exp(r_bh)) ;                           %%/* dL / d b_i */
exp(c_i) + exp(q_h) * (exp(h_i) - exp(h_i(-1))) + ( 1 + exp(r_bh(-1))) * exp(b_i(-1)) / exp(pi) = exp(w_i) *  exp(l_i) + exp(b_i) ;          %%%                                                                                

//main equation!!!
@#if amort 
  @#if amortvar
    rho_amk = 0.5 + rho_am*( exp(y1) / exp(y1(-4) ));
  @#else
    rho_amk = 0.5 + rho_am;
  @#endif
  exp(b_i) = (1-rho_amk)*exp(b_i(-1))/exp(pi) + exp(m_i)*exp(q_h) * (exp(h_i)-(1-delta)*exp(h_i(-1)));      %%%/*BORROWING CONSTRAINT for IH*/
@#else
  (1 + exp(r_bh)) * exp(b_i) = exp(m_i) * exp(q_h(+1)) * exp(h_i) * exp( pi(+1));      %%%/*BORROWING CONSTRAINT for IH*/
@#endif

exp(lambda_e) = (1 - hab_e) / ( exp(c_e) - hab_e * exp(c_e(-1))) ;                                              %%/* dL / d c_e */
exp(lambda_e) * exp(q_k) = exp(s_e) * exp(m_e) * exp (q_k(+1)) * exp(pi(+1)) * (1 - delta) + beta_e * exp(lambda_e(+1)) * ( exp(r_k(+1)) * exp(u(+1)) + exp(q_k(+1))* (1 - delta) 
       - (csi1 * (exp(u(+1)) - 1) + csi2/2 *((exp(u(+1)) - 1)^2))) ;    %%%   /*dL / d k_e*/  
exp(r_k) = csi1 + csi2 *( exp(u) -1 ) ;                                                                         %%/* dL / d u   DEGREE OF UTILIZATION OF CAPITAL*/
exp(r_k) = myalp * exp(a_e) * exp(k_e(-1))^(myalp-1) * exp(u)^(myalp -1) * ( exp(l_ep)^mu * exp(l_ei)^(1-mu) ) ^ (1-myalp) / exp(x) ;  %%%/* MARGINAL PRODUCT OF CAPITAL */
exp(y_e) = exp(a_e) * ( exp(k_e(-1)) * exp(u))^(myalp) * ( exp(l_ep)^mu * exp(l_ei)^(1-mu) ) ^ (1-myalp) ;      %%%/* PRODUCTION FUNCTION*/
exp(w_p) = (1 - myalp) * exp(y_e) * mu / (exp(x) * exp(l_ep) );               %%%/* dL / d l_ep */
exp(w_i) = (1 - myalp) * exp(y_e) * (1 - mu) / (exp(x) * exp(l_ei));          %%%/* dL / d l_ei */
exp(lambda_e) - exp(s_e) *(1 + exp(r_be))  =  beta_e * exp(lambda_e(+1)) * ( 1 + exp(r_be)) / exp(pi(+1)) ;      %%%/* dL / d b_e */
exp(c_e) + (exp(w_p)*exp(l_ep) + exp(w_i)*exp(l_ei)) + ((1 + exp(r_be(-1))) * exp(b_ee(-1)) / exp(pi))  + exp(q_k) * exp(k_e) +
(csi1*(exp(u)  - 1) + csi2/2*(exp(u) -1)^2) * exp(k_e(-1)) = exp(y_e)/exp(x) + exp(b_ee) + exp(q_k) *(1 - delta) * exp(k_e(-1)) ;             %%%      /*BC*/

//main equation!!!
(1 + exp(r_be)) * exp(b_ee) = exp(m_e) * exp( q_k(+1)) * exp( pi(+1)) * (1 - delta) * exp(k_e)  ;       %%%/* BORROWING CONSTRAINT*/

(1 -exp(epsilon_l)) * exp(l_p)+ exp(l_p)^(1 + phi) / exp(w_p) *exp(epsilon_l)/ exp(lambda_p) -
k_w * (exp(pi_wp) - exp(pi(-1))^iota_w * piss^(1 - iota_w)) * exp(pi_wp) + beta_p * exp(lambda_p(+1)) / exp(lambda_p) * k_w * ( exp(pi_wp(+1))  - exp(pi)^( iota_w) * piss^( 1 - iota_w)) * exp(pi_wp(+1))^(2) / exp(pi) = 0;
exp(pi_wp) = exp(w_p) / exp(w_p(-1)) * exp(pi) ;                    %%% /* WAGE INFLATION FOR PATIENT HOUSEHOLDS*/
(1 -exp(epsilon_l)) * exp(l_i) + exp(l_i)^(1 + phi) / exp(w_i) *exp(epsilon_l) / exp(lambda_i) - 
k_w * (exp(pi_wi) - exp(pi(-1))^iota_w * piss^(1 - iota_w)) * exp(pi_wi) + beta_i *  exp(lambda_i(+1)) / exp(lambda_i) * k_w * ( exp(pi_wi(+1))  - exp(pi)^(iota_w) * piss^( 1 - iota_w)) * exp(pi_wi(+1))^(2) / exp(pi) = 0 ;
exp(pi_wi) = exp(w_i) / exp(w_i(-1)) * exp(pi) ;            %%/* WAGE INFLATION FOR IMPATIENT HOUSEHOLDS*/
exp(k) = (1 - delta) * exp(k(-1)) + ( 1 - k_i/2 *(exp(inv) * exp(epsilon_qk)/exp(inv(-1)) - 1 )^2)* exp(inv);            %%%/*AMOUNT OF NEW CAPITAL THAT CGP FIRMS CAN PRODUCE */
exp(q_k) * (1 - k_i/2 * (exp(inv) * exp(epsilon_qk)/ exp(inv(-1)) - 1) ^2 - 
k_i*(exp(inv) * exp(epsilon_qk)/exp(inv(-1)) - 1) * exp(inv) *exp(epsilon_qk)/exp(inv(-1))) + 
beta_e *  exp(lambda_e(+1))/exp(lambda_e) * exp(q_k(+1)) * k_i * ( exp(inv(+1)) * exp(epsilon_qk(+1)) / exp(inv)  - 1 )* exp(epsilon_qk(+1)) * ( exp(inv(+1)) / exp(inv))^2 = 1;  %%% /*REAL PRICE OF CAPITAL*/
1 - exp(epsilon_y) + 
exp(epsilon_y) / exp(x) - 
k_p * (exp(pi) -      (exp(pi(-1))^iota_p * piss^(1 - iota_p) ))* exp(pi) +
beta_p * (exp(lambda_p(+1))/ exp(lambda_p)) *  k_p * (exp(pi(+1)) - (exp(pi)^     iota_p* piss^(1 - iota_p)))* exp(pi(+1)) *(exp(y(+1)) / exp(y)) = 0 ;    %%% /* NON LINEAR PHILLIPS CURVE*/
exp(j_r) = exp(y) *(1 - (1/ exp(x)) - (k_p/2)*(exp(pi) - (exp(pi(-1))^(iota_p) * piss^(1 - iota_p)))^2);           %%%/* RETAILERS PROFITS*/
exp(R_b) = exp(r_t) - k_kb * ( (exp(K_b) / (exp(B)))  - exp(v_b) ) *( exp(K_b) / (exp(B)))^2;     %%%/*FOCS FROM THE MAXIMIZATION PROBLEM OF THE WHOLESALE BRANCH*/
/*exp(r_t) = exp(R_d) ;  *//*NO ARBITRAGE CONDITION WITH THE LENDING FACILITY OF THE CENTRAL BANK*/

@#if vbvar
  @#if y_one_B_zeros
    exp(v_b) =  (vb_ss+0.09) ^( 1-0.92 )  *  ((exp(y1) / exp(y1(-4))-1)*5+1) *  exp( v_b(-1) )^0.92;  
  @#else
    exp(v_b) =  (vb_ss+0.09) ^( 1-0.92 )  *  ((exp(B) / exp(B(-4))-1)*5+1) *  exp( v_b(-1) )^0.92;
  @#endif
@#else
  exp(v_b) =  (vb_ss+0.09) ^( 1-0.92 ) *  exp( v_b(-1) )^0.92;
@#endif

exp(weight_BH) = ( 1 ^( 1-0.94 ) ) * ( exp(y1) / exp(y1(-4)) )^( ( 1 - 0.94 ) * ( - 10 * 1 ) ) * exp( weight_BH(-1) )^0.94; //32
exp(weight_BE) = ( 1 ^( 1-0.92 ) ) * ( exp(y1)/exp(y1(-4)) )^(( 1 - 0.92 ) * ( - 15 * 1  ))* exp( weight_BE(-1) )^0.92; //33

1 - exp(epsilon_bh)/(exp(epsilon_bh) -1)  + exp(epsilon_bh)/(exp(epsilon_bh) -1)  * exp(R_b) / exp(r_bh) - k_bh*(exp(r_bh)/exp(r_bh(-1))  -1 ) * exp(r_bh)/exp(r_bh(-1))  +  
beta_p * (exp(lambda_p(+1)) / exp(lambda_p)) * k_bh * ( exp(r_bh(+1)) / exp(r_bh)  - 1 )  *(( exp(r_bh(+1)) / exp ( r_bh))^2) * (exp(b_h(+1))/ exp(b_h)) = 0 ; 

1 - exp(epsilon_be)/(exp(epsilon_be) -1)   +  exp(epsilon_be)/(exp(epsilon_be) -1)  * exp(R_b) / exp(r_be) - k_be*(exp(r_be)/exp(r_be(-1))  -1 ) * exp(r_be)/exp(r_be(-1))  +  
beta_p * (exp(lambda_p(+1)) / exp(lambda_p)) * k_be * ( exp(r_be(+1)) / exp(r_be)  - 1 )  *(( exp(r_be(+1)) / exp ( r_be))^2) * (exp(b_e(+1))/ exp(b_e)) = 0 ; 

- 1 + exp(epsilon_d)/(exp(epsilon_d) -1)   -  exp(epsilon_d)/(exp(epsilon_d) -1)  * exp(r_t) / exp(r_d) - k_d*(exp(r_d)/exp(r_d(-1))  -1 )  * exp(r_d)/exp(r_d(-1))  +  
beta_p * (exp(lambda_p(+1)) / exp(lambda_p)) * k_d * ( (exp(r_d(+1)) / exp(r_d))  - 1 )   *(( exp(r_d(+1)) / exp ( r_d))^2)   * (exp(d_b(+1))/ exp(d_b)) = 0 ;   

//main equation!!!
exp(j_b)  = + exp(r_bh) * exp(b_h) 
            + exp(r_be) * exp(b_e)
            + exp(r_t)  * exp(B)*eta
            - exp(r_d)  * exp(d_b) 
            - delta_b*exp(K_b)- k_kb/2 * ( ((exp(K_b) / (exp(wB)))  - exp(v_b) ) ^2) * exp(K_b);      % I added delta_b here
              
exp(K_b) * exp(pi) = (1 - delta_b) * exp(K_b(-1))  + (pi_share)*exp(j_b(-1)) ;             %%%          /*BANK CAPITAL ACCUMULATION*/
gamma_b * exp(d_b)  = gamma_p * exp(d_p) ;                %%%                      
gamma_b * exp(b_h)  = gamma_i * exp(b_i) ;                %%%
gamma_b * exp(b_e)  = gamma_e * exp(b_ee);                %%%
exp(B)*(1 - eta)  =  exp(d_b) + exp(K_b)  + eps_K_b;       %%%/*BALANCE SHEET CONSTRAINT  B = D + K_b*/
(1 + exp(r_t)) = (1 + r_t_ss)^(1 - phi_R) * (1 + exp(r_t(-1)))^(phi_R) *
                 ((exp(pi)/piss)^(phi_pi) * (exp(y1)/exp(y1(-1)))^phi_y)^(1 - phi_R) *(1+e_r_t) ;    %%%

//exp(y)  = exp(c) + 
//          /*exp(q_k)* */    ( exp(k) - (1 - delta) *exp(k(-1))) ;
//          /*+ exp(k_e(-1)) * (csi1*(exp(u) -1 ) + (csi2/2)*(exp(u) -1)^(2)) + */
//          /*delta_b* exp(K_b(-1)) / exp(pi)  */
         

//          /*- (k_p/2) * (exp(pi) /   - exp(pi(-1))^iota_p * piss^(1-iota_p)) ^(2)*exp(y) */
//          /*- k_d/2  * ( (exp(r_d)/exp(r_d(-1))-1)^2)   * exp(r_d) *exp(d_b)  */
//          /*- k_be/2 * ( (exp(r_be)/exp(r_be(-1))-1)^2) * exp(r_be)*exp(b_e)  */
//          /*- k_bh/2 * ( (exp(r_bh)/exp(r_bh(-1))-1)^2) * exp(r_bh)*exp(b_h);  */ 
exp(c_pr)= gamma_p * exp(c_p) + gamma_i * exp(c_i) + gamma_e * exp(c_e);
exp(c_g) = pcg2pcp*exp(c_p);
exp(c)   = gamma_p * exp(c_p) + gamma_i * exp(c_i) + gamma_e * exp(c_e) + exp(c_g)  ;%%
h        = gamma_p * exp(h_p) + gamma_i * exp(h_i) ;           %%
exp(k)   = gamma_e * exp(k_e) ;
exp(B)*(1 - eta)   = (exp(BH) + exp(BE)) ;

@#if weightvar
  wbhk = (1+wbh);
@#else
  wbhk = (1+wbh) * (1+5*(exp(y1) / exp(y1(-4))-1));
@#endif

exp(wB)*(1 - eta)   = ((1+wbhk)*exp(BH) + exp(BE));
exp(BH)  = gamma_b * exp(b_h);
exp(BE)  = gamma_b * exp(b_e);
exp(D)   = gamma_p * exp(d_p);
exp(y)   = gamma_e * exp(y_e) ;
exp(J_b) = gamma_b * exp(j_b);  
gamma_e * exp(l_ep) = gamma_p * exp(l_p);
gamma_e * exp(l_ei) = gamma_i * exp(l_i);

exp(tax) = exp(c_g) + (eta)*( exp(d_b(-1)-exp(d_b) - exp(K_b)) + exp(K_b(-1)) )*(exp(r_t(-1))) ; //gov.equation
exp(spr_b)    = 0.5*exp(r_bh) + 0.5*exp(r_be) - exp(r_d); //Int.rate spread

exp(PIW)            = ( exp(w_p) + exp(w_i) )  / ( exp(w_p(-1)) + exp(w_i(-1)) ) * exp(pi);
exp(y1)             = exp(c) +    1     * (exp(k)-(1-delta)*exp(k(-1))) ;
                   //13

/***************************************************************+ EXOGENOUS PROCESS*****************************************************************************/

exp(epsilon_z)     = 1 - rho_e_z   *    1               + rho_e_z   * exp(epsilon_z(-1))    + e_z;        %%
exp(a_e)           = 1 - rho_e_ae  *    1               + rho_e_ae  * exp(a_e(-1))          + nega;       %%
exp(epsilon_h)     = 1 - rho_e_h   *    1               + rho_e_h   * exp(epsilon_h(-1))    + e_h;        %%

@#if ltvvar
  @#if y_one_B_zeros
    exp(m_i)           = (1-rho_mi)    *  (0.75 + m_i_ss + myltv)  *( 1 - (exp(y1) / exp(y1(-4))-1)*5) + rho_mi    * exp(m_i(-1))          + e_mi;       %% /*y1 was B*/
    exp(m_e)           = (1-rho_me)    *  (m_e_ss + myltv) *( 1 - (exp(y1) / exp(y1(-4))-1)*5)           + rho_me    * exp(m_e(-1))          + e_me;       %%
  @#else
    exp(m_i)           = (1-rho_mi)    *  (0.75 + m_i_ss + myltv)  *( 1 - (exp(B) / exp(B(-4))-1)*5) + rho_mi    * exp(m_i(-1))          + e_mi;       %% /*y1 was B*/
    exp(m_e)           = (1-rho_me)    *  (m_e_ss + myltv) *( 1 - (exp(B) / exp(B(-4))-1)*5)           + rho_me    * exp(m_e(-1))          + e_me;       %%
  @#endif
@#else
  exp(m_i)           = (1-rho_mi)    *  (0.75 + m_i_ss + myltv)  + rho_mi    * exp(m_i(-1))          + e_mi;       %%
  exp(m_e)           = (1-rho_me)    *  (m_e_ss + myltv)            + rho_me    * exp(m_e(-1))          + e_me;       %%
@#endif

//exp(m_e)           = (1-rho_me)    *  m_e_ss            + rho_me    * exp(m_e(-1))          + e_me;       %%
exp(epsilon_d)     = (1-rho_e_d)   * epsilon_d_ss       + rho_e_d   * exp(epsilon_d(-1))    + e_d;        %%
exp(epsilon_be)    = (1-rho_e_be)  * epsilon_be_ss      + rho_e_be  * exp(epsilon_be(-1))   + e_be;       %%
exp(epsilon_bh)    = (1-rho_e_bh)  * epsilon_bh_ss      + rho_e_bh  * exp(epsilon_bh(-1))   + e_bh;       %%
exp(epsilon_qk)    =  1-rho_e_qk   *    1               + rho_e_qk  * exp(epsilon_qk(-1))   + e_qk;       %%
exp(epsilon_y)     = (1-rho_e_y)   * epsilon_y_ss       + rho_e_y   * exp(epsilon_y(-1))    + e_y;        %%
exp(epsilon_l)     = (1-rho_e_l)   * epsilon_l_ss       + rho_e_l   * exp(epsilon_l(-1))    + e_l;        %%
exp(eps_K_b)       = (1-rho_eps_K_b)*    1              + rho_eps_K_b* exp(eps_K_b(-1))     + e_eps_K_b;  %%

y_obs = y;
c_obs = c;
i_obs = inv;
q_h_obs = q_h;
r_t_obs = r_t;
r_be_obs = r_be;
r_bh_obs = r_bh;
r_d_obs = r_d;
b_e_obs = b_e;
b_h_obs = b_h;
d_b_obs = d_b;
pi_obs = pi;
w_i_obs = w_i;
kb_rat = exp(K_b)/exp(B); 

end;