[pth] = fileparts(which('vwi'));
home_dir = char(textread([pth '\home_dir.txt'],'%s'));

%% Prompt for Directory to process and define all subdirectories.
uiwait(msgbox('Please select directory you would like to analyze.','VWI'));
proc_dir = uigetdir(home_dir, 'Select directory you would like to analyze..');
dir_proc = dir(proc_dir);

size_dir = size(dir_proc,2);

for ii = 1:1:size_dir-2
    subname = dir_proce.name(ii);
    subdir = [proc_dir '\' subname];
    excel_name = dir([subdir,'\' subname '_statistics.xlsx']);
    data = xlsread(excel_name);
    
