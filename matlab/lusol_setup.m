function lusol_setup(install)
% LUSOL_SETUP setup script for lusol
%
% Usage:
%   >> cd [lusol directory]
%   >> lusol_setup('install') % to create load_lusol file in userpath
%   >> lusol_setup            % to add lusol directly to path
%
% After lusol_setup('install'), the following command will add lusol to
% Matlab's path:
%   >> load_lusol
%

if nargin == 1 && strcmp(install,'install')
  % get userpath directory
  mypath = strsplit(userpath,':');
  mypath = mypath{1};
  % create directory if it does not exist
  % TODO: make this better
  [~,~,~] = mkdir(mypath)
  % open file
  fid = fopen([mypath '/load_lusol.m'],'w');
  % write file
  fprintf(fid,'function load_lusol\n');
  fprintf(fid,'  addpath(''%s'')\n',pwd);
  fprintf(fid,'end\n');
  % close file
  fclose(fid);
  %keyboard
else
  addpath(pwd);
end

end
