% clear;clc

%myvars = {'v_b','u', 'r_k', 'r_d', 'r_t', 'r_be', 'r_bh', 'spr_b', 'c_p', 'c_i', 'c_e' 'c', 'inv', 'y', 'k_e', 'q_h', 'q_k', 'h_i', 'd_p', 'b_h', 'b_e', 'l_i', 'l_p', 'j_r', 'K_b'};
myvars = {'y', 'c', 'inv', 'r_d', 'r_t', 'r_be', 'r_bh', 'spr_b', 'd_p', 'b_h', 'b_e', 'l_i', 'l_p', 'h_i', 'kb_rat', 'K_b'};
myvars_long = {'Output', 'Consumption', 'Investment', 'Deposit rate', 'Policy rate', 'Loan rate (E)', ...
               'Loan rate (HH)', 'Interest rate spread', 'Deposit', 'Loan (HH)', 'Loan (E)', 'Labor (imp.)', ...
               'Labor (p.)', 'Housing (HH)', 'CAR', 'Bank capital'};
% myvars = {'r_t' 'r_be' 'r_bh' 'l_i' 'l_p' 'c' 'inv' 'q_h' 'h_i' 'h_p' 'd_p' 'spr_b' };
% shock = {'e_ae','e_r_t'};
% shock = {'e_ae', 'e_r_t', 'e_h','e_mi', 'e_me','e_l', 'e_qk', 'e_y', 'e_bh', 'e_be', 'e_d', 'e_eps_K_b'};

packagemFiles({'define1.mod','main_body.mod','initval.mod','shock_later.mod'});
dynare combined nostrict
for i=1:numel(myvars)
    eval(['myfirstmod.' myvars{i} '=' myvars{i} ';']);
end

packagemFiles({'define1.mod','main_body.mod','initval.mod','shock_early.mod'});
dynare combined nostrict
for i=1:numel(myvars)
    eval(['mysecondmod.' myvars{i} '=' myvars{i} ';']);
end

t=1:1:numel(myfirstmod.([myvars{1}]));


%   figure;
%   for i=1:numel(myvars)
%     
%     subplot(5,5,i)
%     plot(t,myvarmod.([myvars{i} '_' shock{1}]),'b', t, myvarmod.([myvars{i} '_' shock{2}]),'k');
%     ax=gca;
%     ax.GridColor = [0 0 0];
%     ax.GridLineStyle = '--';
%     title(myvars{i});
%     ytickformat('%.2f');
%     grid;
%     xlim([0 numel(t)]);
%   end
%   title_ = ['IRF to ' shock{j}];
%   suptitle(title_)
%   leg=legend({'Var', 'Fixed'});
%   newPosition = [0.5 0.05 0 0];
%   newUnits = 'normalized';
%   set(leg,'Units', newUnits,'Position', newPosition, 'NumColumns', 2);
  


% for j=1:numel(shock)
  figure;
  for i=1:numel(myvars)
    
    subplot(4,4,i)
    plot(t,myfirstmod.([myvars{i}]),'b--', t, mysecondmod.([myvars{i}]),'b');
    ax=gca;
    ax.GridColor = [0 0 0];
    ax.GridLineStyle = '--';
    title(myvars_long{i});
    ytickformat('%.2f');
    grid;
    xlim([0 numel(t)]);
  end
% title(detvar);
%   suptitle(title_)
  leg=legend({'Late (4 quarter ahead)', 'Early (8 quarter ahead)'});
  newPosition = [0.5 0.05 0 0];
  newUnits = 'normalized';
  set(leg,'Units', newUnits,'Position', newPosition, 'NumColumns', 2);
  fpath = 'C:\Users\Sugar2\Documents\dyn\Master\figuress';
  filename = [detvar '_det'];
  saveas(gcf, fullfile(fpath, filename), 'png');
% end

