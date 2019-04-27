% myvars = {'kb_rat','u', 'r_k', 'r_d', 'r_t', 'r_be', 'r_bh', 'spr_b', 'c_p', 'c_i', 'c_e' 'c', 'inv', 'y', 'k_e', 'q_h', 'q_k', 'h_i', 'd_p', 'b_h', 'b_e', 'l_i', 'l_p', 'j_r', 'K_b'};
myvars = {'y', 'c', 'inv', 'r_d', 'r_t', 'r_be', 'r_bh', 'spr_b', 'd_p', 'b_h', 'b_e', 'l_i', 'l_p', 'h_i', 'kb_rat', 'K_b'};
myvars_long = {'Output', 'Consumption', 'Investment', 'Deposit rate', 'Policy rate', 'Loan rate (E)', ...
               'Loan rate (HH)', 'Interest rate spread', 'Deposit', 'Loan (HH)', 'Loan (E)', 'Labor (imp.)', ...
               'Labor (p.)', 'Housing (HH)', 'CAR', 'Bank capital'};
shock = {'e_ae','e_r_t'};
mypar_grid = struct('myltv', [-0.05 0 0.15], ...
                    'k_kb', [0 5 15]);
mypar_grid_names = fieldnames(mypar_grid);

for k = 1:numel(mypar_grid_names)
irfs=[];
mypar = mypar_grid_names(k);
myparss = [mypar{1} 's'];
eval([myparss '= mypar_grid.(mypar_grid_names{k});']);
mypars =eval(myparss);

first_time = 1;
for i=1:length(mypars)
    if first_time
%         set_param_value(mypar{1},mypars(i));
          dynare combined noclearall nostrict;
        irfs=[irfs oo_.irfs];
        first_time = 0;
    else
        set_param_value(mypar{1},mypars(i));
        info = stoch_simul(var_list_);
          irfs=[irfs oo_.irfs];
        if info
          disp(['Computation fails for delta = ' num2str(mypar{1})]);
        end
    end
end

t=1:1:numel(irfs(1).a_e_e_ae);

figure;

for j=1:numel(shock)
  for i=1:numel(myvars)
    subplot(4,4,i)
    p = plot(t, irfs(1).([myvars{i} '_' shock{j}]),'b--',t,irfs(2).([myvars{i} '_' shock{j}]),'b',t,irfs(3).([myvars{i} '_' shock{j}]),'ro-','MarkerSize',2);
    set(p, {'linewidth'}, {1});
    title([([myvars_long{i}])]); grid; 
    ax=gca;
    ax.GridColor = [0 0 0];
    ax.GridLineStyle = '--';
    if i == 5
      leg=legend({[mypar{1} ' = ' num2str(mypars(1))], [mypar{1} ' = ' num2str(mypars(2))], [mypar{1} ' = ' num2str(mypars(3))]}, 'Interpreter', 'none');
      newPosition = [0.5 0.05 0 0];
      newUnits = 'normalized';
      set(leg,'Units', newUnits,'Position', newPosition, 'NumColumns', 3);
    end
  end
fpath = 'C:\Users\Sugar2\Documents\dyn\Master\figuress';
filename = [mypar{1} '_'  shock{j} '_' ind_ltvvar];
saveas(gcf, fullfile(fpath, filename), 'png');  
end



end
% suptitle([mypar{1} ' grid'])
% leg=legend({[mypar{1} ' = ' num2str(mypars(1))], [mypar{1} ' = ' num2str(mypars(2))], [mypar{1} ' = ' num2str(mypars(3))]}, 'Location', 'SouthOutside', 'boxoff');
% newPosition = [0.5 0.2 0 0];
% newUnits = 'normalized';
% set(leg,'Units', newUnits, 'NumColumns', 3);
  
%end