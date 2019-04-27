

initval;
vb_ss = 0;
//m_e_ss = 0;
m_i_ss = 0;
rho_am = 0;
wbh=0;
pi = log(piss);
pi_wp = pi;
pi_wi = pi;
PIW = pi;

  @#if amort
  rho_amk = 0.3;
  @#endif

% Shocks
q_k = 0; % eq.25
epsilon_z = log(1); % eq.59
a_e = log(1); % eq.60
epsilon_h = log(1); % eq.61
m_i = log(0.75); % eq.62
m_e = log(m_e_ss); % eq.63
epsilon_d = log(epsilon_d_ss); % eq.64
epsilon_bh = log(epsilon_bh_ss); % eq.65
epsilon_be = log(epsilon_be_ss); % eq.66
epsilon_qk = log(1); % eq.67
epsilon_y = log(epsilon_y_ss); % eq.68
epsilon_l = log(epsilon_l_ss); % eq.69  
eps_K_b = log(1); % eq.70

% Interest rates
r_t = log(r_t_ss); % eq.41
R_b = r_t; % - k_kb*0.03*1; % eq.28 Kb/B != v_b
r_bh = R_b + epsilon_bh; % eq.32
r_be = R_b + epsilon_be; % eq.33
r_d = r_t + epsilon_d; % eq.34
spr_b = log(0.5*exp(r_bh)+0.5*exp(r_be)-exp(r_d));
% Utilization
u = log(((((exp(q_k) - (1-beta_e*(1+exp(r_be)))/(1+exp(r_be))*exp(m_e+q_k)*(1-delta)-beta_e*exp(q_k)*(1-delta))/beta_e-csi1+csi2/2)*2/csi2))^(0.5)); 
r_k = log(csi1 + csi2 *( exp(u)-1));

% 
x = log(exp(epsilon_y)/(exp(epsilon_y)-1));

% l_ei
l_ei =-0.08 + log(((1+exp(r_bh+epsilon_h)*J/(1+exp(r_bh))/(1-beta_i-exp(pi+m_i)*(1-beta_i*(1+exp(r_bh))/exp(pi))/(1+exp(r_bh))))*(-(1-exp(epsilon_l))/(exp(epsilon_l))))^(1/(mu+phi)));
//l_ei = log(((1+exp(r_bh+epsilon_h)*J/(1+exp(r_bh))/(1-beta_i-exp(pi+m_i)*(1-beta_i*(1+exp(r_bh))/exp(pi))/(1+exp(r_bh))))*(-(1-exp(epsilon_l))/(exp(epsilon_l))))^(1/(mu+phi)));
l_i = l_ei;
% l_ep
%dbh = (exp(r_bh) + exp(r_t)*eta/(1-eta) - (1+delta/pi_share))/(exp(r_d)-1 - delta/pi_share);
%dbh = (exp(r_bh) + exp(r_t)*eta/(1-eta) - (delta+delta/pi_share+0.0081*k_kb))/(exp(r_d)-delta - delta/pi_share-0.0081*k_kb);
%dbe = (exp(r_be) + exp(r_t)*eta/(1-eta) - (1+delta/pi_share))/(exp(r_d)-1 - delta/pi_share);
%dbe = (exp(r_be) + exp(r_t)*eta/(1-eta) - (delta+delta/pi_share+0.0081*k_kb))/(exp(r_d)-delta - delta/pi_share-0.0081*k_kb);

%m1 = (exp(r_d)-delta*(1-pi_share)/pi_share)*(exp(r_bh) + exp(r_t)*eta/(1-eta) - (1+delta/pi_share))/(exp(r_d)-1 - delta/pi_share) + delta*(1-pi_share)/pi_share;
%m2 = (exp(r_d)-delta*(1-pi_share)/pi_share)*(exp(r_bh) + exp(r_t)*eta/(1-eta) - (1+delta/pi_share))/(exp(r_d)-1 - delta/pi_share) + delta*(1-pi_share)/pi_share;
%m3 = (1-1/exp(x));
%m0 = mu*(1-myalp)/exp(x);
%mm1 = J*exp(epsilon_h+m_i)/(1-beta_i-exp(m_i)*(1-beta_i*(1+exp(r_bh))/exp(pi)/(1+exp(r_bh)))); % s_i/lambda_i
%mm2 =  myalp*(1-delta)*exp(m_e-r_k-u-x)/(1+exp(r_be));

%mmm1 = (1-mu)*(1-myalp)/exp(x)/(1+exp(r_bh)/(1+exp(r_bh))*J*exp(epsilon_h+m_i)/(1-beta_i-exp(m_i)*(1-beta_i*(1+exp(r_bh))/exp(pi)/(1+exp(r_bh)))));

%ms = (mu*(1-myalp)/exp(x) + J*exp(epsilon_h+m_i)/(1-beta_i-exp(m_i)*(1-beta_i*(1+exp(r_bh))/exp(pi)/(1+exp(r_bh))))*J*exp(epsilon_h+m_i)/(1-beta_i-exp(m_i)*(1-beta_i*(1+exp(r_bh))/exp(pi)/(1+exp(r_bh))))*(1-mu)*(1-myalp)/exp(x)/(1+exp(r_bh)/(1+exp(r_bh))*J*exp(epsilon_h+m_i)/(1-beta_i-exp(m_i)*(1-beta_i*(1+exp(r_bh))/exp(pi)/(1+exp(r_bh)))))+(exp(r_d)-delta*(1-pi_share)/pi_share)*(exp(r_bh) + exp(r_t)*eta/(1-eta) - (1+delta/pi_share))/(exp(r_d)-1 - delta/pi_share) + delta*(1-pi_share)/pi_share*myalp*(1-delta)*exp(m_e-r_k-u-x)/(1+exp(r_be))+(1-1/exp(x)));

//l_ep=log((-(1-exp(epsilon_l))*mu*(1-myalp)/(exp(epsilon_l)*(mu*(1-myalp)/exp(x) + J*exp(epsilon_h+m_i)/(1-beta_i-exp(m_i)*(1-beta_i*(1+exp(r_bh))/exp(pi)/(1+exp(r_bh))))*J*exp(epsilon_h+m_i)/(1-beta_i-exp(m_i)*(1-beta_i*(1+exp(r_bh))/exp(pi)/(1+exp(r_bh))))*(1-mu)*(1-myalp)/exp(x)/(1+exp(r_bh)/(1+exp(r_bh))*J*exp(epsilon_h+m_i)/(1-beta_i-exp(m_i)*(1-beta_i*(1+exp(r_bh))/exp(pi)/(1+exp(r_bh)))))+(exp(r_d)-delta*(1-pi_share)/pi_share)*(exp(r_bh) + exp(r_t)*eta/(1-eta) - (delta+delta/pi_share+0.0081*k_kb))/(exp(r_d)-delta - delta/pi_share-0.0081*k_kb) + delta*(1-pi_share)/pi_share*myalp*(1-delta)*exp(m_e-r_k-u-x)/(1+exp(r_be))+(1-1/exp(x)))*exp(x)))^(1/(1+phi)));
l_ep=log((-(1-exp(epsilon_l))*mu*(1-myalp)/(exp(epsilon_l)*(mu*(1-myalp)/exp(x) + J*exp(epsilon_h+m_i)/(1-beta_i-exp(m_i)*(1-beta_i*(1+exp(r_bh))/exp(pi)/(1+exp(r_bh))))*J*exp(epsilon_h+m_i)/(1-beta_i-exp(m_i)*(1-beta_i*(1+exp(r_bh))/exp(pi)/(1+exp(r_bh))))*(1-mu)*(1-myalp)/exp(x)/(1+exp(r_bh)/(1+exp(r_bh))*J*exp(epsilon_h+m_i)/(1-beta_i-exp(m_i)*(1-beta_i*(1+exp(r_bh))/exp(pi)/(1+exp(r_bh)))))+(exp(r_d)-delta*(1-pi_share)/pi_share)*(exp(r_bh) + exp(r_t)*eta/(1-eta) - (1+delta/pi_share))/(exp(r_d)-1 - delta/pi_share) + delta*(1-pi_share)/pi_share*myalp*(1-delta)*exp(m_e-r_k-u-x)/(1+exp(r_be))+(1-1/exp(x)))*exp(x)))^(1/(1+phi)));
l_p = l_ep;
% k

k_e = log(((exp(r_k+x)/(myalp*exp(a_e)))^(1/(myalp-1))*exp(l_ep)^mu*exp(l_ei)^(1-mu))/exp(u));
k = k_e;

% y_e y
y_e = log(exp(a_e)*exp(k_e+u)^myalp*(exp(l_ep)^mu * exp(l_ei)^(1-mu))^(1-myalp));
y = y_e;

% wages
w_i = log((1-mu)*(1-myalp)*exp(y-x)/exp(l_ei));
w_p = log(mu*(1-myalp)*exp(y-x)/exp(l_ep));

% c
c_p = log(-exp(w_p)*(1-exp(epsilon_l))/(exp(epsilon_l)*exp(l_ep)^phi));
c_i = log(-exp(w_i)*(1-exp(epsilon_l))/(exp(epsilon_l)*exp(l_ei)^phi));

b_ee = log(exp(m_e) * exp(q_k) * exp(pi) * (1 - delta) * exp(k_e)/(1 + exp(r_be)));
b_e = b_ee;
BE = b_e;

c_e = log(exp(y_e)/exp(x) + exp(b_ee) + exp(q_k) *(1 - delta) * exp(k_e) - ((exp(w_p)*exp(l_ep) + exp(w_i)*exp(l_ei)) + ((1 + exp(r_be)) * exp(b_ee) / exp(pi))  + exp(q_k) * exp(k_e) + (csi1*(exp(u)  - 1) + csi2/2*(exp(u) -1)^2) * exp(k_e)));  

% lambda
lambda_i = - c_i;
lambda_p = - c_p;
lambda_e = - c_e;



s_i = log((exp(lambda_i)-beta_i*exp(lambda_i)*(1+exp(r_bh)))/(1+exp(r_bh)))+0.05;
s_e = log((exp(lambda_e) - beta_e*exp(lambda_e)*(1+exp(r_be)))/(1+exp(r_be)));

//h_p = exp(lambda_p)*(1-beta_p)/(exp(lambda_i)*(1-beta_i-exp(pi+s_i+m_i)))/(1+exp(lambda_p)*(1-beta_p)/(exp(lambda_i)*(1-beta_i-exp(pi+s_i+m_i))));
h_p = exp(lambda_p)*(1-beta_p)/(exp(lambda_i)*(1-beta_i-exp(pi+s_i+m_i)))/(1+exp(lambda_p)*(1-beta_p)/(exp(lambda_i)*(1-beta_i-exp(pi+s_i+m_i))))-0.016;
h_i = log(1- exp(h_p))-0.45;
//h_i = log(1- exp(h_p));

q_h = -log(exp(lambda_p+h_p)*(1-beta_p)/(J*exp(epsilon_h)));


@#if amort
  b_i = log(exp(m_i)*exp(q_h) * (delta*exp(h_i))/rho_amk);      %%%/*BORROWING CONSTRAINT for IH*/
@#else
  b_i = log(exp(pi+h_i+q_h+m_i)/(1+exp(r_bh)));
@#endif

b_h = b_i;
BH= b_h;

% d = dbh*bh + dbe*be
d_p = log((exp(r_bh) + exp(r_t)*eta/(1-eta) - (1+delta/pi_share))/(exp(r_d)-1 - delta/pi_share)*exp(b_h) + (exp(r_be) + exp(r_t)*eta/(1-eta) - (1+delta/pi_share))/(exp(r_d)-1 - delta/pi_share)*exp(b_e));
d_b = d_p;
D = d_b;

B = log((exp(b_h)+exp(b_e))/(1-eta));
wB = log(((1+wbh)*exp(b_h)+exp(b_e))/(1-eta));
K_b = log(exp(B)*(1-eta)-exp(d_p)-eps_K_b);

v_b = log(0.09);
j_r = log(exp(y) *(1 - (1/ exp(x))));
j_b  = log(exp(r_bh) * exp(b_h)+ exp(r_be) * exp(b_e)+ exp(r_t)  * exp(B)*eta- exp(r_d)  * exp(d_b)- exp(K_b)- k_kb/2 * ( ((exp(K_b) / (exp(wB)))  - exp(v_b) ) ^2) * exp(K_b));  
J_b = j_b;


inv = log(delta*exp(k));

weight_BH = log(1);
weight_BE = log(1);

c_pr= log(gamma_p * exp(c_p) + gamma_i * exp(c_i) + gamma_e * exp(c_e));
c_g = log(pcg2pcp*exp(c_p));
tax = c_g;
c = log(gamma_p * exp(c_p) + gamma_i * exp(c_i) + gamma_e * exp(c_e) + exp(c_g));
y1 = log(exp(c) + delta*exp(k));
end;

check;
steady;