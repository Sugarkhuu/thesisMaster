
val = '0.1';

fname_early = 'shock_later.mod';
fileID_early = fopen(fname_early,'w');

fprintf(fileID_early,'\n');
fprintf(fileID_early,'shocks;\n');
fprintf(fileID_early,['var ' detvar ';\n']);
fprintf(fileID_early,['periods 5:12;\n']);
fprintf(fileID_early,['values ' val ';\n']);
fprintf(fileID_early,'end;\n');
fprintf(fileID_early,'simul(periods=24);\n');
  
fclose(fileID_early);


fname_late = 'shock_early.mod';
fileID_late = fopen(fname_late,'w');

fprintf(fileID_late,'\n');
fprintf(fileID_late,'shocks;\n');
fprintf(fileID_late,['var ' detvar ';\n']);
fprintf(fileID_late,['periods 9:16;\n']);
fprintf(fileID_late,['values ' val ';\n']);
fprintf(fileID_late,'end;\n');
fprintf(fileID_late,'simul(periods=24);\n');
  
fclose(fileID_late);