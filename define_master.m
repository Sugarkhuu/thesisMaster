function define_master(fname, varargin)

if isempty(varargin)
  define_struct = struct('weightvar', 0, ...
                       'ltvvar', 0,    ...
                       'vbvar', 0,     ...
                       'kbfree', 0,    ...
                       'amort', 0,     ...
                       'amortvar', 0,  ...
                       'kldfree', 0,   ...
                       'freebank', 0);
else
  define_struct = varargin{:};
end

atdef = '@#define';
fn = fieldnames(define_struct);                                      
fileID = fopen(fname,'w');
  for i = 1:numel(fn)
    
    fprintf(fileID,[atdef ' ' '%s = %i\n'], fn{i}, define_struct.(fn{i}));
  
  end
fclose(fileID);


end