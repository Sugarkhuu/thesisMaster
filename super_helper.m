%% Fixed

weightvars = [0 0 0 0 0 0 0 0 0 0];
ltvvars    = [0 1 0 0 0 0 0 0 0 0];
vbvars     = [0 0 0 0 0 0 0 0 0 0];
kbfrees    = [0 0 0 0 0 0 0 0 0 0];
amorts     = [0 0 0 0 0 0 0 0 0 0];
amortvars  = [0 0 0 0 0 0 0 0 0 0];
kldfrees   = [0 0 0 0 0 0 0 0 0 0];
freebanks  = [0 0 0 0 0 0 0 0 0 0];
y_one_B_zeros  = [0 1 0 0 0 0 0 0 0 0];

for i = 1:2
def1 = i;
defineka = struct('weightvar', weightvars(def1), ...
                       'ltvvar', ltvvars(def1),    ...
                       'vbvar', vbvars(def1),     ...
                       'kbfree', kbfrees(def1),    ...
                       'amort', amorts(def1),     ...
                       'amortvar', amortvars(def1),  ...
                       'kldfree', kldfrees(def1),   ...
                       'freebank', freebanks(def1), ...
                       'y_one_B_zeros', y_one_B_zeros(def1));
                
nameka = 'define1.mod';
define_master(nameka, defineka);
packagemFiles({'define1.mod','main_body.mod','initval.mod','sto_simul.mod'});

if ltvvars(def1) == 1
  ind_ltvvar = 'ltvvar';
else 
  ind_ltvvar = '';
end

param_grid_combined;

end

%% Fixed vs Variable

weightvars = [0 0 0 0 0 0 0 0 0 0];
ltvvars    = [0 0 0 1 0 0 1 1 0 0];
vbvars     = [0 1 0 0 1 1 0 0 0 0];
kbfrees    = [0 0 0 0 0 0 0 0 0 0];
amorts     = [0 0 0 0 0 0 0 0 0 0];
amortvars  = [0 0 0 0 0 0 0 0 0 0];
kldfrees   = [0 0 0 0 0 0 0 0 0 0];
freebanks  = [0 0 0 0 0 0 0 0 0 0];
y_one_B_zeros  = [0 0 0 0 1 0 1 0 0 0];

for i = 1:((numel(weightvars)-2)/2)
  
nameka = 'define1.mod';
def1 = 2*i;
defineka = struct('weightvar', weightvars(def1), ...
                       'ltvvar', ltvvars(def1),    ...
                       'vbvar', vbvars(def1),     ...
                       'kbfree', kbfrees(def1),    ...
                       'amort', amorts(def1),     ...
                       'amortvar', amortvars(def1),  ...
                       'kldfree', kldfrees(def1),   ...
                       'freebank', freebanks(def1), ...
                       'y_one_B_zeros', y_one_B_zeros(def1));
define_master(nameka, defineka);

nameka = 'define2.mod';
def2 = 2*i-1;
defineka = struct('weightvar', weightvars(def2), ...
                       'ltvvar', ltvvars(def2),    ...
                       'vbvar', vbvars(def2),     ...
                       'kbfree', kbfrees(def2),    ...
                       'amort', amorts(def2),     ...
                       'amortvar', amortvars(def2),  ...
                       'kldfree', kldfrees(def2),   ...
                       'freebank', freebanks(def2), ...
                       'y_one_B_zeros', y_one_B_zeros(def2));
define_master(nameka, defineka);

if ltvvars(def1) == 1
  ind_ltvvar = 'ltvvar';
else 
  ind_ltvvar = '';
end
if vbvars(def1) == 1
  ind_vbvar = 'vbvar';
else 
  ind_vbvar = '';
end

if y_one_B_zeros(def2) == 1
  ind_y_one_B_zero = 'y_one_B_zero';
else 
  ind_y_one_B_zero = '';
end

fix_var;

end

%% Deterministic 

weightvars = [0 0];
ltvvars    = [1 1];
vbvars     = [1 1];
kbfrees    = [0 0];
amorts     = [0 0];
amortvars  = [0 0];
kldfrees   = [0 0];
freebanks  = [0 0];
y_one_B_zeros  = [0 0];

for i = 1:numel(weightvars)

nameka = 'define1.mod';
def2 = i;
defineka = struct('weightvar', weightvars(def2), ...
                       'ltvvar', ltvvars(def2),    ...
                       'vbvar', vbvars(def2),     ...
                       'kbfree', kbfrees(def2),    ...
                       'amort', amorts(def2),     ...
                       'amortvar', amortvars(def2),  ...
                       'kldfree', kldfrees(def2),   ...
                       'freebank', freebanks(def2), ...
                       'y_one_B_zeros', y_one_B_zeros(def2));
define_master(nameka, defineka);

detvars = {'m_i_ss','vb_ss'};
detvar = detvars{i};
gen_shocks;

my_det_sim;

end

%% Properties

weightvars = [0 0 0 0 0 0 0 0 0 0];
ltvvars    = [0 0 0 1 0 0 1 1 0 0];
vbvars     = [0 1 0 0 1 1 0 0 0 0];
kbfrees    = [0 0 0 0 0 0 0 0 0 0];
amorts     = [0 0 0 0 0 0 0 0 0 0];
amortvars  = [0 0 0 0 0 0 0 0 0 0];
kldfrees   = [0 0 0 0 0 0 0 0 0 0];
freebanks  = [0 0 0 0 0 0 0 0 0 0];
y_one_B_zeros  = [0 0 0 0 1 0 1 0 0 0];

for i = 1
  
nameka = 'define1.mod';
def1 = i;
defineka = struct('weightvar', weightvars(def1), ...
                       'ltvvar', ltvvars(def1),    ...
                       'vbvar', vbvars(def1),     ...
                       'kbfree', kbfrees(def1),    ...
                       'amort', amorts(def1),     ...
                       'amortvar', amortvars(def1),  ...
                       'kldfree', kldfrees(def1),   ...
                       'freebank', freebanks(def1), ...
                       'y_one_B_zeros', y_one_B_zeros(def1));
define_master(nameka, defineka);

property;

end



