clear;
clc;

repoDir = fileparts(mfilename('fullpath'));
if isempty(repoDir)
    repoDir = pwd;
end
addpath(genpath(repoDir));

main_intervalrobust;
