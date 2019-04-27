function packagemFiles(FileNames)

% varargin(1:2:end) = lower(varargin(1:2:end));
%   cellIx = cellfun(@iscell,varargin);
%   varargin(cellIx) = cellfun(@(x){{x}},varargin(cellIx));
%   opts = struct(varargin{:});
  
  
    N=length(FileNames);
    if N<1
        return;
    else
        TragetFile = ['combined.mod'];
        [hTargetFile,ErrMsg] = fopen(TragetFile,'w');
    end
    for i=1:N
        [hSourceFile,ErrMsg] = fopen(FileNames{i},'r');
        if hSourceFile < 0
            disp([varargin{i},'is not found!']);
            fclose(hSourceFile);
            return; 
        end
        Data = fread(hSourceFile);
        fclose(hSourceFile);
        fwrite(hTargetFile, Data);
%         fprintf(hTargetFile,'\n');
    end
    fclose(hTargetFile);
end