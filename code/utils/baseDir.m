function base_dir = baseDir()
%BASEDIR Returns the base directory used in the project which is defined as
% $HOME/Documents/phd/ 
%   
%   Inputs:
%
%   Return:
%    - base_dir, base directory for the project, $HOME/Documents/phd/.
base_dir = join([getenv("HOME"), "Documents/phd"], '/');
end