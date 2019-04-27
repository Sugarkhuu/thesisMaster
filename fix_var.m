% clear;clc

%myvars = {'v_b','r_k', 'r_d', 'r_t', 'r_be', 'r_bh', 'spr_b', 'c_p', 'c_i', 'c_e' 'c', 'inv', 'y', 'k_e', 'q_h', 'q_k', 'h_i', 'd_p', 'b_h', 'b_e', 'l_i', 'l_p', 'j_r', 'K_b'};
myvars = {'y', 'c', 'inv', 'r_d', 'r_t', 'r_be', 'r_bh', 'spr_b', 'd_p', 'b_h', 'b_e', 'l_i', 'l_p', 'h_i', 'kb_rat', 'K_b'};
myvars_long = {'Output', 'Consumption', 'Investment', 'Deposit rate', 'Policy rate', 'Loan rate (E)', ...
               'Loan rate (HH)', 'Interest rate spread', 'Deposit', 'Loan (HH)', 'Loan (E)', 'Labor (imp.)', ...
               'Labor (p.)', 'Housing (HH)', 'CAR', 'Bank capital'};
shock = {'e_ae','e_r_t'};

packagemFiles({'define1.mod','main_body.mod','initval.mod','sto_simul.mod'});
dynare combined nostrict
for j=1:numel(shock)
  for i=1:numel(myvars)
    eval(['myvarmod.' myvars{i} '_' shock{j} '=' myvars{i} '_' shock{j} ';']);
  end
end
packagemFiles({'define2.mod','main_body.mod','initval.mod','sto_simul.mod'});
dynare combined nostrict

for j=1:numel(shock)
  for i=1:numel(myvars)
    eval(['myfixmod.' myvars{i} '_' shock{j} '=' myvars{i} '_' shock{j} ';']);
  end
end

t=1:1:numel(myvarmod.([myvars{i} '_' shock{1}]));

for j=1:numel(shock)
  figure;
  for i=1:numel(myvars)
    
    subplot(4,4,i)
    p = plot(t,myvarmod.([myvars{i} '_' shock{j}]),'b--', t, myfixmod.([myvars{i} '_' shock{j}]),'b');
    set(p, {'linewidth'}, {1.5});
    ax=gca;
    ax.GridColor = [0 0 0];
    ax.GridLineStyle = '--';
    title(myvars_long{i});
    ytickformat('%.2f');
    grid;
    xlim([0 numel(t)]);
    if i ==1
      if isempty(ind_y_one_B_zero)
        leg=legend({'Variable rule', 'Fixed rule'});
      else
        leg=legend({'Output', 'Loan'});
      end 
      newPosition = [0.5 0.05 0 0];
      newUnits = 'normalized';
      set(leg,'Units', newUnits,'Position', newPosition, 'NumColumns', 2);
    end 
  end
%   title_ = ['IRF to ' shock{j}];
%   suptitle(title_)
  fpath = 'C:\Users\Sugar2\Documents\dyn\Master\figuress';
  filename = [ind_ltvvar  ind_vbvar '_' ind_y_one_B_zero '_' shock{j}];
  saveas(gcf, fullfile(fpath, filename), 'png');
end




