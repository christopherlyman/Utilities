function [final_table] = DTIgrad(Astruct)

% ==========================================
% HEADER
% ==========================================

% Program Name: DTIgrad.m

% Author:   Jonathan A.D. Farrell ( JonFarrell@jhu.edu )
%           F.M. Kirby Research Center for Functional Brain Imaging
%           Kennedy Krieger Institute
%           Baltimore, MD 21218

% Creation Date:  August 2, 2005
% Last Modified:  December 29, 2007

% HOW TO ACKNOWLEDGE / CITE: 		
% http://godzilla.kennedykrieger.org/~jfarrell/software_web/

% LICENSE AND DISCLAIMER INFORMATION:  	
% http://godzilla.kennedykrieger.org/~jfarrell/gpl.txt

disp('****************************************************************************************');
disp('DTIgrad.m, figures out the DTI gradient table');
disp('for the Philips MRI scanner, based on information about your scan');
disp('Copyright (C) 2005  Jonathan A.D. Farrell');
disp(' ');
disp('The DTIgrad.m software was kindly provided by ');
disp('Jonathan A.D. Farrell (Johns Hopkins University, Balitmore, Maryland USA)')
disp('The development of this software was supported by NIH/NCRR grant RR15241 ');
disp('to the Kennedy Krieger Institute.');
disp('For information on how to acknowledge / cite please visit: ');
disp('http://godzilla.kennedykrieger.org/~jfarrell/software_web/ ');
disp(' ');
disp('This program is free software; you can redistribute it and/or modify');
disp('it under the terms of the GNU General Public License as published by')
disp('the Free Software Foundation; either version 2 of the License, or');
disp('(at your option) any later version.');
disp('For license information please visit')
disp('http://godzilla.kennedykrieger.org/~jfarrell/gpl.txt');
disp(' ');
disp('This program is NOT intended for clinical use. USE AT YOUR OWN RISK !');
disp(' ');
disp('This program is distributed in the hope that it will be useful,');
disp('but WITHOUT ANY WARRANTY; without even the implied warranty of');
disp('MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the');
disp('GNU General Public License for more details.');
disp('You should have received a copy of the GNU General Public License along');
disp('with this program; if not, write to the Free Software Foundation, Inc.,');
disp('51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.');
disp('****************************************************************************************');
% ========================================
% PURPOSE:  
% ========================================

% This program is meant for data acquired using the F.M. Kirby center 
% Philips MR scanners. The program adjusts the gradient table (from the 
% source code) to account for several imaging options (like gradient 
% overplus, slice orientation, slice angulation, foldover direction, 
% patient orientation,Philips software release, and software upgrades to 
% the Kirby center scanners etc) that affect the orrientation of the 
% diffusion encoding gradients on the scanner.  The code also adjusts the 
% gradient tables for the rotational component of  coregistration if motion
% correction was performed with AIR or FSL FLIRT.  Both of these affect 
% the gradient tables that are ultimately used to calculate the diffusion 
% tensor 

% When using DTIstudio for fiber tracking, you DO NOT need to flip any 
% eigenvectors.  

% If you acquired data on a differient Philips MR scanenr (NOT at the F.M. 
% Kirby Research Center)please use DTI_gradient_table_creator_Achieva_RelX.m

% =====================
% History of Updates:
% =====================

% August 3, 2005 | Added a few more display comments to the angulation
% correction section

% August 4, 2005 | Validated the results for several test datasets used
% throughout the rotation_ovp.m and DTI_gradient_transform.m validation
% process. Also added a check statement to see if user has UserWrite
% privilages to the path that the .grad file will be written. (by default,
% this is the path of the .par / .PAR file entered). Error message is
% displayed if UserWrite = 0.

% August 4, 2005 | removed typo end statement from line 163. 

% September 22, 2005 | added a line to make sure we don't normalize the
% mean DWI volume direction if it exists in the Angulation corrected table.
% By convention this is labled with [100,100,100] as the direction

% September 28, 2005 | Modified to account for case when specified par file
% has no path, ie Matlab function is called from location of par file.
% Modified to accept string '3' input for scanner in addition to '3.0' and
% 3 as a numeric.

% September 29, 2005 | Added an additional error statement for foldover
% direction that does not match fat_shift options.

% October 10, 2005 | Fixed bug in 3.0 T invert-Jones30 code section.
% and cleaned up some definitions for the variable that hold the
% information about the VMS or Windows OS usage.  I also provided a matlab
% based .air file readin function (scanair.m).  This replaces the call to
% the scanair function on a unix system.  This change should make
% DTIgrad.m operating system independent (UNIX vs
% windows)

% February 3, 2006 | Added FH as one of the accepted options for fold over
% direction% (ie. preperation direction in the V4 par files).  Note that 
% the FH foldover direction is not compatible with the SENSE headcoil.

% March 11, 2006 | Fixed a bug in the parse Header section.  It would crash
% when a lower case v or upper case V were in the patient name. Craig Jones
% fixed code in readrec, I updated my code.

% July 29, 2006 | You can now specify usere-defined as a grad_choice. You 
% must also provide an ASCII text file that has the gradient directions.
% These directions should be EXACTLY as entered on the Philips scanner.
% Changed from function inputs to a structured array input format.  Major
% changes to validation of inputs was required.

% July 31, 2006 | Added a variable to the structured array to allow the
% path to the transformation files from motion correction (i.e., the
% xfmpath) to be different from the path to the par file (i.e., the path)

% August 2, 2006 | added fullfile and strtrim changes to 
% registration_correction function this should clear up the problem of 
% white space at the end of file names and simplify the code a little

% August 9, 2006 | Incorporated Bennett's changes to the checklist 
% validation step

% August 30, 2006 | Removed xfmpath from structure array input validation.  
% If xfmpath is not a field, or is a field and has en empty array value 
% A.xfmpath = [], then I search for the transformation files in the same 
% path as the par file.  If you specify A.xfmpath = '/C/bla/bla/bla', then
% I search in this directory for the transformation files.  Also, if the 
% A.supplied_grad_file is not specified in the structured array, I set it 
% equal to an empty array. Made correction to registration_correction 
% section to make sure the xfm_air_directory has a slash on the end of the 
% filename. If it does not, I add one, if it does,then I don't add anything.

% September 13, 2006 | Added detailed information about the licence 
% agreement and disclaimer.

% September 18, 2006 | found bug in statement elseif 
% findstr(position_flag,'Left') found bug, had firndstr. 
% Have fixed it to read findstr

% November 2, 2006 | I nolonger append the 0,0,0 vector to the begining of 
% the user-defined table. I trust that the user will provide it in the 
% Astruct.supplied_grad_file text file

% December 13, 2006 | I added Rel_2.0 and Rel_2.1 to the list of supported
% Philips software releases.  I also added the option to sort or not to 
% sort the images.  if you sort the images then the b=0 volume as the 
% first volume and the mean DWI volume as the last.  If you choose not to
% sort the images, the b=0 is the 2nd last, then mean DWI is the last.
% The new input paramter is Astruct.sort_images = 'y', or 'n'. Since
% sorting the images is customary, I will set a default of 'y' if you fail
% to provide an input for Astruct.sort_images.The Kirby center scanners
% will not be on the same software version.  The 3T will move to R2.1 while
% the 1.5 T will remain at R11.x for some time.

% December 14, 2006 | To avoid confusion I changed the variable 'path' to
% 'parpath'.  The path variable in Matlab can have special meaning.  I also
% added some support for CATNAP. The name of the .air or .xfm files may be
% different from the name of the par file, so long as the name is passed in
% the variable Astruct.xfmname.  If this entry is empty or missing, I
% default to the par name.

% December 19, 2006 | Added a validation section to make sure the slice
% orientation and foldover are legal.  for each slice orientation there are
% only 2 legal foldover options (out of the 3 foldover possibilities of AP,
% FH and RL). 

% March 14, 2007 | The 3T scanner in the Kirby center has been  upgraded
% to Rel_2_1.  I have updated the code accordingly.  At this time, the 1.5T
% scanner will remain as Rel_11

% April 4, 2007 | If a required input field in Astruct is not defined, 
% I tell you what is missing.  This is done in the input validation section
% I also replaced the isfield statement which used an array input.  This
% was done to work with some versions of Matlab 7.

% April 11, 2007 | I am discovered a bug in how we were interpreting the
% sort images button on the scanner software.  The history is as follows.
% for V3 and V4 par files, if the sort images button was checked, then the
% REC file is "volume" ordered.  This means that slices 1 through N for the
% 1st volume are written, then slices 1 through N for the 2nd volume, etc.
% Additionally, the b=0 volume is listed as the first volume.  I do not
% have examples of what happened when the sort images button was not
% checked for V3.  Colleagues suggest that this option did not work correctly.  It
% will not be supported for V3 par files. DO NOT SUPPLY THE 0,0,0 (ie. b=0 )
% volume in your user supplied gradient table.  I will append it as necessary.

% >>>> A MAJOR CHANGE occurred with the upgrade to V4.1 PAR files.
% In V4.1 par files, if the sort images button is checked on the scanner,
% the REC file is "volume" ordered. This means that slices 1 through N for the
% 1st volume are written, then slices 1 through N for the 2nd volume, etc.
% HOWEVER, the b=0 volume is now listed as the 2nd last volume (after the mean
% DWI volume), or as the last volume (if the mean DWI volume is not
% present).  In V4.1 par files, if the sort images button is NOT checked,
% then the data is written out as "slice ordered".  This means that slices 1 
% for volumes 1 through M are written, followed by slice 2 for volumes 1
% through M. Additionally, the b=0 volume is listed as the first sub volume.

% June 15, 2007 | I made a correction to the input validation for V3 par
% files.  In the section that checks for the definition of foldover,
% patient_position and patient_orientation I had a misplaced line which
% caused a failure if the par file was V3.  I also improved the error
% message that is displayed if any of the required inputs are not provided.

% SUMMARY:

% Version   sort_images     order        b=0 position in REC file
%   V3          yes           volume       first volume
%   V3          no           NOT SUPPORTED
%   V4          yes           volume       first volume
%   V4          no            slice        first sub volume for slices
%   V4.1        yes           volume       last or second last volume
%   V4.1        no            slice        first sub volume for slices

% IMPACT on gradient tables for DTIStudio:

% This code writes out gradient tables for DTIStudio.  The gradient table
% provided will reflect the position of the b=0 volume in the REC file (it
% will be first, last or second last depending on the table above).  For
% the case of V4.1 sort_images = y, the b=0 volume is last.  For the case
% of V4.1 sort_images = 'n', the b=0 volume must be listed first to work
% with DTIStudio and image sequence "slice by slice" must be selected 

% July 20, 2007 | I corrected a small bug with the rotation matrices for
% slice angulation.  I had multiplied 3 matrices in the incorrect order. 
% A colleague (Harsh Agarwal) pointed this out while aligning different
% MRI contrasts using the angulation parameters and the transformation 
% matrices given in the Philips document.  
%I originally had Tang = Tfh*Tap*Trl
%    which is now Tang = Trl*Tap*Tfh;
%I originally had rev_Tang = rev_Trl*rev_Tap*rev_Tfh; 
%which is now     rev_Tang = rev_Tfh*rev_Tap*rev_Trl;
% I double checked the Philips code and this seems to be correct now.
% I also double checked the impact on fiber tracking. The fiber tracking
% looks good in both instances (even though the gradient tables are
% slightly different).  If 2 angulation values are zero (i.e. [AP,FH,RL] =
% [0,0,20], then the old and new equations give the same result.  Only if
% two or more elements are non zero is the result different.  I did some
% testing with very large angulations of 20 degrees [20,20,0], [20,0,20] 
% and [0,20,20]and found that the fiber tracking results were almost 
% indistinguishable. THIS FIX ONLY affects yes overplus and 
% user-defined gradient tables. No overplus (and therefor Jones30) tables 
% are not subject to slice angulation changes.

% July 31, 2007 | Added option of entering gradient tables using the
% user-file option on the Kirby patch.  Note, this is completely different
% to using the user-defined option.  The user-file gradient tables play out
% as no overplus tables, wheras the user-defined tables are like
% yes-overplus tables.  If you select grad_choice = 'user-file' then  you
% need to provide a gradient table in a text file.

% November 13, 2007 | Added support for the new yes overplus high gradient table
% for release 2.5. Philips changed one direction in this table to remove a 
% duplicate direction that was in previous versions.  Added par file version 4.2
% support and updated the parsing operations for reading in the par file

% December 29, 2007 | Fixed bug regarding Releaes 2.5 yes overplus gradient tables
% I had forgotten to add Rel_2_5 to the list of checks when I assign the space as
% XYZ (pre release 1.2) or LPH (release 1.2 and later, including release 2.5).
% The release 2.5 yes overplus tables were incorrectly assigned as XYZ space,
% when they should in fact be LPH space.  Incorrect colormaps will have their red and 
% green colors interchanged and fiber tracking will be incorrect. 



% ==========================================
% GENERAL COMMENTS ABOUT INPUTS:
% ============================================

% If you are using V4 or more recent par files, I can get the foldover, 
% patient,orientation, patient_position and all slice orientation and 
% slice angulation parameters from the .par and .PAR files.  

% However, if  you are using V3 par files, then you MUST enter foldover, 
% patient_orientation and patient_position as none of these
% parameters are in the .par or .PAR file.  I do however get the slice 
% orientation and slice angulation parameters from the .par and .PAR files.

% ==============================
% Example of REQUIRED input
% ================================
%A.par_file = '/g2/asif/data/LO/year2/060626_pm/060626_pm_14_1.par'; 
% input the path and filename for your par file.  example.
% /g2/data/sample.par 

%A.didREG = 'n'; 
%'y' or 'n' If you did coregistration (AKA motion correction' enter 'y', 
% The the .air (AIR transformation files) or .xfm (FLIRT transformation files) 
% may be in the /g2/data/ path (i.e., the same path as the par file) or may
% be in a seperate xfmpath (specified in the structured input).  If you
% fail to provide a seperate xfmpath, or you define it as an empty array,
% then the default is xfmpath = path (i.e., I look in the path of the par
% file). If you did not do coregistration of your data, then put 'n' for 
% didREG and you don't need the .air or .xfm files
% For example.  in the /g2/data/ directory will would see
% sample.par
% sample_001.air
% sample_002.air
% ....
% sample_00N.air where N is the last volume of your acquisition.

% If didREG is 'y', then what I do is search for any files that are sample*.air or sample*xfm
% assuming that you did EITHER AIR or FLIRT registration.  I then sort them into
% alphabetical order (001 then 002 then 003 .... 00N), then use the
% transformation in each .air or .xfm file to correct the corresponding
% line in the diffusion gradient table.

%A.writeGRAD = 'n'; 
% 'y' or 'n' will write results to file called <par_file>.grad

%A.grad_choice = 'jones30'; 

% 'yes-ovp-low', 'yes-ovp-medium','yes-ovp-high',''no-ovp-low',
% 'no-ovp-medium','no-ovp-high','jones30','invert-jones30','user-defined', 
% or 'user-file'.
% This specifies the gradient table used during the DTI acquisition.  NOTE: NOT ALL OF
% THESE OPTIONS ARE AVAILABLE, DEPENDING ON THE DATE OF YOUR SCAN. The
% program is built to catch any conflicts between known dates of upgrades,
% the information in your par file and your input for grad_choice.
% NOTE.  If you specify grad_choice = 'user-defined', or 'user-file' then 
% you must also provide an ASCII text file that has the gradient directions listed
% EXACTLY as entered on the Philips scanner, do not include the 0,0,0
% direction as this will be added automatically.

%A.scanner = '1.5'; 
%'1.5' or 1.5 or '3.0', '3' or 3.0 VERY important to get correct tables on
% KIRBY scanners

%A.fat_shift = 'P'; 
% 'A','P','R','L','H','F'. This is the fat shift direction

%A.sort_images = 'y'; % 'n'

% =======================================================
% ADDITIONAL REQUIRED inputs if par version is V3
% ======================================================
%A.foldover = 'AP'; 
% 'AP','RL' or 'FH' for foldover.  This is the phase encoding 
% (AKA preperation direction)

%A.patient_orientation = 'sp'; 
% 'sp','pr','rd' or 'ld' 
% This is for supine,prone,right decubitus and left decubitus.

%A.patient_position = 'hf'; % 'hf' or 'ff' This is for head first or feet first

% ==========================================================
% ADDITIONAL REQUIRED input if grad_choice = 'user-defined'
% ==========================================================
%A.supplied_grad_file = []; 
% provide an ASCII text file that has the gradient directions listed
% EXACTLY as entered on the Philips scanner. do not include the 0,0,0
% direction as this will be added automatically.

% ===================================
% Inputs specific to the Kirby center
% ===================================
%Astruct.flag = 'KIRBY'; 
% 'KIRBY' or 'OTHER'
%Astruct.release = 'Rel_11.x'; % or 'Rel_2.1' or 'Rel_2.5'
% the exact release depends on the software upgrade date on the 1.5T or 3T
% scanners

% ====================================
% ====================================
% START OF INPUT VALIDATION
% ====================================
% ====================================

q = {'par_file','didREG','writeGRAD','grad_choice','scanner','fat_shift'};
for j=1:length(q)
    checklist(j) = isfield(Astruct,q{j});
end

for ii=1:length(q);
    checklist = isfield(Astruct,q{ii});
    if (checklist == 0)
        disp('You failed to provide the REQUIRED inputs')
        disp('One or more of par_file, didREG, writeGRAD, grad_choice, scanner, or fat_shift are not defined in Astruct')
        error([q{ii} ' is not defined']);
    end
end

par_file = Astruct.par_file;
if exist(par_file) ~= 2
    disp('You need to enter a valid par file')
    error(['Sorry...The par file : ' par_file ' you specified does not exist'])
end

[parpath,name,ext] = fileparts(par_file);
if (~isfield(Astruct,'xfmname'))
% If empty, defualt to the par name
elseif isempty(Astruct.xfmname)
% If not a field , default to the par name    
elseif (isfield(Astruct,'xfmname') & ~isempty(Astruct.xfmname))
    % overide the name with the xfmname
    name = Astruct.xfmname;
end

% get some information from the .par / .PAR file
[nrows, ncols, nslices, nechoes, ndynamics, nphases, A, header, par_version] = int_getPARinfo(par_file);
par_version = num2str(par_version)

if ~(strcmpi(par_version,'3') | strcmpi(par_version,'4') | strcmpi(par_version,'4.1') | strcmpi(par_version,'4.2'))
    error('unsupported and unknown par file vesion, must be 3, 4, 4.1 or 4.2')
end
% Count the number of volumes reported in the par file.  This number will
% be useful for distinguishing between varieties of scans
nvolumes = size(A,1)/nslices;  

% ======================================
% Special V3 Par file version validation step
% =======================================
if (strcmpi(par_version,'4') | strcmpi(par_version,'4.1') | strcmpi(par_version,'4.2'))
    
elseif (strcmpi(par_version,'3'))
    test = {'foldover','patient_orientation','patient_position'};
    for ii=1:length(test);
        checklist = isfield(Astruct,test{ii});
        if (checklist == 0)
            disp('Since your par file is older than V4, (it is most likely a V3), you must enter')
            disp('inputs for the foldover, patient_orientation and patient_position in Astruct')
            error([test{ii} ' is not defined']);
        end
    end
end

% =====================================
% Continue with input validation steps
% =====================================
didREG = Astruct.didREG;
if (strcmpi(didREG,'y') | strcmpi(didREG,'n'))
else
    disp('Did you perform coregistration of your DTI data?')
    error('Please enter the didREG option as y or n')
end

writeGRAD = Astruct.writeGRAD;
if (strcmpi(writeGRAD,'y') | strcmpi(writeGRAD,'n'))
else
    disp('Do you want to write out the output of this program as a .grad file?')
    disp('The .grad file will contain the gradient table and is in a format ')
    disp('that is used in DTIstudio and dtiproc based analysis')
    error('Please enter the didREG option as y or n')
end 

grad_choice = Astruct.grad_choice;
if (strcmpi(grad_choice,'yes-ovp-high') | strcmpi(grad_choice,'yes-ovp-medium') | strcmpi(grad_choice,'yes-ovp-low')...
    |strcmpi(grad_choice,'no-ovp-high') | strcmpi(grad_choice,'no-ovp-medium') | strcmpi(grad_choice,'no-ovp-low')...
    |strcmpi(grad_choice,'jones30') | strcmpi(grad_choice,'invert-jones30') | strcmpi(grad_choice,'user-defined')...
    |strcmpi(grad_choice,'user-file'))
else
    error('Please enter: yes-ovp-high,yes-ovp-medium,yes-ovp-low,no-ovp-high,no-ovp-medium,no-ovp-low,jones30,invert-jones30,user-defined or user-file for first input ')
end  

if (strcmpi(grad_choice,'user-defined') | strcmpi(grad_choice,'user-file'))
    if isfield(Astruct,'supplied_grad_file')
        supplied_grad_file = Astruct.supplied_grad_file;
        if exist(supplied_grad_file) ~= 2
            disp('You need to enter a valid gradient table file')
            error(['Sorry...The gradient table file : ' supplied_grad_file ' you specified does not exist'])
        end
    else
        error('Use of the user-defined option or user-file option requires that you provide a file with the directions you entered at the scanner')
    end
else
    supplied_grad_file = [];
end

fat_shift = Astruct.fat_shift;
if (strcmpi(fat_shift,'A') | strcmpi(fat_shift,'P') | strcmpi(fat_shift,'R')...
    |strcmpi(fat_shift,'L') | strcmpi(fat_shift,'F') | strcmpi(fat_shift,'H'))
else
    error('Please enter the fat shift direction as A, P, R, L, F,or H ')
end

scanner = Astruct.scanner;
if ~isnumeric(scanner)
    if (strcmpi(scanner,'1.5') | strcmpi(scanner,'3.0') | strcmpi(scanner,'3'))
        if strcmpi(scanner,'3')
            scanner = '3.0';
        end
    else
        error('Please enter the scanner choice as 1.5, 3 or 3.0 as strings or numerics')
    end
else
    if (scanner == 1.5)
        scanner = '1.5';
    elseif (scanner == 3.0)
        scanner = '3.0';
    else 
        error('Please enter the scanner choice as 1.5 or 3.0 as strings or numerics')
    end
end

% If not a field , default to the par path
if (~isfield(Astruct,'xfmpath'))
    xfmpath = parpath;
% If empty, defualt to the par path
elseif isempty(Astruct.xfmpath)
    xfmpath = parpath;
elseif (isfield(Astruct,'xfmpath') & ~isempty(Astruct.xfmpath))
    if(exist(Astruct.xfmpath) == 7)
        xfmpath = Astruct.xfmpath;
    else
        error(['the xfmpath you provided: ' Astruct.xfmpath ' is not a valid directory'])
    end
end

% If not a field, default to 'y'
if (~isfield(Astruct,'sort_images'))
    sort_images = 'y';
elseif isempty(Astruct.sort_images)
    sort_images = 'y';
elseif (isfield(Astruct,'sort_images') & ~isempty(Astruct.sort_images))
    sort_images = Astruct.sort_images;
    if (strcmpi(sort_images,'y') | strcmpi(sort_images,'n'))
    else
        disp('Do you want to sort the images?')
        error('Please enter the sort_images option as y or n')
    end
end

% ============================================
% Gather inputs from V4 and more recent par files
% ============================================
if (strcmpi(par_version,'4') | strcmpi(par_version,'4.1') | strcmpi(par_version,'4.2'))
    
    % get the patient position
    position_flag = header.patient_position;
    foldover_flag = header.preparation_direction;
    if findstr(position_flag,'Head')
        patient_position = 'hf';
    elseif findstr(position_flag,'Feat')
        patient_position = 'ff';
    end
    
    % get the patient_orientation
    if findstr(position_flag,'Supine')
        patient_orientation = 'sp';
    elseif findstr(position_flag,'Prone')
        patient_orientation = 'pr';
    elseif findstr(position_flag,'Right')
        patient_orientation = 'rd';
    elseif findstr(position_flag,'Left') 
        patient_orientation = 'ld';
    end
    
    % get the slice orientation
    slice_flag = unique(A(:,26))
    if size(slice_flag) == 1
        if slice_flag == 1
            slice_orientation = 'tra';
        elseif slice_flag == 2
            slice_orientation = 'sag';
        elseif slice_flag == 3
            slice_orientation = 'cor';
        end
    else
        error('You have more than one slice orientation listed in the par / PAR file')
    end
    
    % get the foldover direction
    if findstr(header.preparation_direction,'Anterior-Posterior')
        foldover = 'AP';
    elseif findstr(header.preparation_direction,'Right-Left')
        foldover = 'RL';
    elseif findstr(header.preparation_direction,'Feet-Head')
        foldover = 'FH';
    end     
elseif (strcmpi(par_version,'3')) % If older than a V4 (i.e., V3) 
    % get the foldover direction
    foldover = Astruct.foldover;
    if (strcmpi(foldover,'AP') | strcmpi(foldover,'RL') | strcmpi(foldover,'FH'))
    else
        error('Please enter the foldover direction as AP, RL, or FH')
    end

    % get the patient orientation
    patient_orientation = Astruct.patient_orientation;
    if (strcmpi(patient_orientation,'sp') | strcmpi(patient_orientation,'pr') | strcmpi(patient_orientation,'rd')...
        |strcmpi(patient_orientation,'ld'))
    else
        error('Please enter the patient orientation as sp,pr,rd or ld')
    end

    % get the patient position
    patient_position = Astruct.patient_position;
    if (strcmpi(patient_position,'hf') | strcmpi(patient_position,'ff'))
    else
        error('Please enter the patient position as hf or ff')
    end  
    
    % get the slice orientation for V3 par files
    slice_flag = unique(A(:,20));
    if size(slice_flag) == 1
        if slice_flag == 1
            slice_orientation = 'tra';
        elseif slice_flag == 2
            slice_orientation = 'sag';
        elseif slice_flag == 3
            slice_orientation = 'cor';
        end
    else
        error('You have more than one slice orientation listed in the par / PAR file')
    end
    
end

% ======================================
% Check that foldover and fat_shift agree
% =========================================
if (strcmpi(foldover,'AP'))
    if ~(strcmpi(fat_shift,'A') | strcmpi(fat_shift,'P'))
        error(['The fat shift ' fat_shift ' and foldover ' foldover ' you provided are not compatible'])
    end
elseif (strcmpi(foldover,'RL'))
    if ~(strcmpi(fat_shift,'R') | strcmpi(fat_shift,'L'))
        error(['The fat shift ' fat_shift ' and foldover ' foldover ' you provided are not compatible'])
    end
elseif (strcmpi(foldover,'FH'))
    if ~(strcmpi(fat_shift,'F') | strcmpi(fat_shift,'H'))
        error(['The fat shift ' fat_shift ' and foldover ' foldover ' you provided are not compatible'])
    end
end

% =========================================================
% Check that the foldover and slice orientation are legal
% =========================================================

if strcmpi(slice_orientation,'tra')
    if ~(strcmpi(foldover,'AP') | strcmpi(foldover,'RL'))
        error('The foldover and slice orientation are not compatible')
        disp('With transverse slices you can choose AP or RL foldover')
    end
     
elseif strcmpi(slice_orientation,'cor')
    if ~(strcmpi(foldover,'FH') | strcmpi(foldover,'RL'))
        error('The foldover and slice orientation are not compatible')
        disp('With coronal slices you can choose RL or FH foldover')
    end
        
elseif strcmpi(slice_orientation,'sag')
    if ~(strcmpi(foldover,'FH') | strcmpi(foldover,'AP'))
        error('The foldover and slice orientation are not compatible')
        disp('With sagital slices you can choose AP or FH foldover')
    end
end

ap = header.angulation(1);
fh = header.angulation(2);
rl = header.angulation(3);
scan_date = header.examination_date;

% =========================================
% ==========================================
% END INPUT PARAMETER VALIDATION
% =========================================
% =========================================

% ==========================================
% CALL TO ANGULATION CORRECTION
% =========================================
DTI_studio = 'y';
[input_table,Ang_corrected_table] = angulation_correction(grad_choice,patient_orientation,patient_position,slice_orientation,foldover,fat_shift,ap,fh,rl,'NA.txt','n',DTI_studio,scan_date,scanner,nvolumes,supplied_grad_file,sort_images,par_version); 
clear input_table;

if strcmpi(didREG,'y')
	% ======================================================
	% CALL TO REGISTRATION CORRECTION IF NEEDED
	% ====================================================
	[final_table] = registration_correction(xfmpath, name, Ang_corrected_table); 
else
    final_table = Ang_corrected_table;
end

if strcmpi(writeGRAD,'y')
    % ==========================
    % WRITE RESULTS TO FILE  
    % ==========================
	filename = [name '.grad'];
	% write the out the table to a .grad file if needed
	i = 0:1:size(final_table,1)-1;
	t = [i;permute(final_table,[2 1])];
    
    % check to see if you have user write privilages
    [stat,mess] = fileattrib(xfmpath);
    if mess.UserWrite == 1
		fid = fopen(fullfile(xfmpath,filename),'w');
		fprintf(fid,['%2.0f:  %5.4f, %5.4f, %5.4f\n'],t);
		fclose(fid);
        disp('Calculations are complete')
    else
        disp(['You do not have UserWrite privilages to the destination ' xfmpath ' so cannot write the .grad file'])
        error(['Please give yourself UserWrite privilages to ' xfmpath ])
    end
else
    disp('Calculations are complete')
end






% THE FUNCTIONS BELOW ARE CALLED BY THE CODE ABOVE
% DO NOT MAKE CHANGES BELOW....!!!

% ============================================================
% ============================================================
% START OF ANGULATION CORRECTION
% ============================================================
% ============================================================
function [input_table,Ang_corrected_table] = angulation_correction(grad_choice,patient_orientation,patient_position,slice_orientation,foldover,fat_shift,ap,fh,rl,filename,doWrite,DTI_studio,scan_date,scanner,nvolumes,supplied_grad_file,sort_images,par_version); 

% ==============================
% HEADER
% ==============================

% Name:     angulation_correction.m (2nd generation of rotation_ovp.m)

% Author:   Jonathan Farrell, 
%           F.M. Kirby Research Center for Functional Brain Imaging
%           Room G25
%           Kennedy Krieger Institute
%           Baltimore, MD 21218

% Email:         JonFarrell@jhu.edu

% Creation Date: February 2, 2005

% ================================
% UPDATES:
% ================================

% June 16, 2005   | Cleaned up some code to reduce the number of variables used
% , by eliminatint the ovp variable and incoorporating the information into
% the grad_choice variable. Added functionality and comments dealing with the 1.5T
% upgrade and jones30 and invert-jones30 implementation on the 1.5 and 3T.

% July 2, 2005.  | Input for scanner can now be numeric 1.5 or 3.0 or string
% 1.5 or string 3.0. Makes sense to make it this way.  Cleaned up some
% display statements. 

% July 6, 2005.  | Discovered that the table used on the VMS 1.5 T system by Joe Gillen had
% an incorrect gradient direction.  The offending direction is number 17.
% The first element in the vector had its sign reversed, by accident.
% I made the changes to rotation_ovp so that the correct table to analyze
% the DATA from before the upgrade will be generated.  In summary, the table used before the
% upgrade was meant to be the Jones30, but had one sign transposed by
% accident.  I now generate the correct Jones30 -1 table which has the
% correct direction for the DWI data before the upgrade.

% July 7, 2005.  | Moved the scan date display line above the point where I
% change the meaning of scan_date.

% August 1, 2005 | MAJOR REVISIONS -- Updated rotation_ovp to read in information from the par
% files when possible. This will make rotation_ovp easier to use.  Revised the order of input parameters to be in the
% most essential first, then others if required based on par version.  NEED
% TO FIX COMMENTS.  THIS VERSION IS UNTESTED.  The m file is now called
% angulation_correction.m

% August 3, 2005.  | Did the testing using DTIgrad.  see
% notes above.

% July 29 2006 | added supplied_grad_file

% November 2, 2006 | In the angulation_correction section I was normalizing all vectors, 
% but if the 0,0,0 vector was provided (when using the user-defined opeion), then a NaN was produced. 
% I now skip over any 0,0,0 or 100,100,100 vectors during the normalization step.  I also no
% longer append the 0,0,0 vector to the begining of the user-defined table.  I trust that the
% user will provide it.

% December 13, 2006 | I added Rel_2.0 and Rel_2.1 to the list of supported
% Philips software releases.  I also added the option to sort or not to 
% sort the images.  if you sort the images then the b=0 volume as the 
% first volume and the mean DWI volume as the last.  If you choose not to
% sort the images, the b=0 is the 2nd last, then mean DWI is the last.
% The new input paramter is Astruct.sort_images = 'y', or 'n'. Since
% sorting the images is customary, I will set a default of 'y' if you fail
% to provide an input for Astruct.sort_images. The Kirby center scanners
% will not be on the same software version.  The 3T will move to R2.1 while
% the 1.5 T will remain at R11.x for some time.

% March 14, 2007 | The 3T scanner in the Kirby center has been  upgraded
% to Rel_2_1.  I have updated the code accordingly.  At this time, the 1.5T
% scanner will remain as Rel_11

% April 11, 2007 | do not enter the b=0 volume as part of your user defined
% table

% July 20, 2007 | corrected Tang and rev_Tang (see comment for main
% program)

% November 13, 2007 | added support for Rel 2.5

% December 29, 2007 | Fixed bug regarding Releaes 2.5 yes overplus gradient tables
% I had forgotten to add Rel_2_5 to the list of checks when I assign the space as
% XYZ (pre release 1.2) or LPH (release 1.2 and later, including release 2.5).
% The release 2.5 yes overplus tables were incorrectly assigned as XYZ space,
% when they should in fact be LPH space.  Incorrect colormaps will have their red and 
% green colors interchanged and fiber tracking will be incorrect. 

% define the dates when the scanners switched from VMS to Windows Operating
% systems.  This is involved in the implementation of the Jones 30 gradient
% table.  In the VMS era, the table has an error.
VMS_to_windows_upgrade_15_date = datenum([2005,04,26]);
VMS_to_windows_upgrade_30_date = datenum([2003,11,24]);

% ============ waiting for upgrade to Rel 2.1 =========================
% Define the dates when the scanners switched from Release 11 to Release
% 2.1.  This is important for the overplus = yes gradient tables.  The
% coordinate space of these tables is known to change from Release 1.2 to
% Release 1.5.  
%Rel_11_to_Rel_2_1_upgrade_15_date = datenum([2006,12,25]); 
Rel_11_to_Rel_2_1_upgrade_30_date = datenum([2007,03,14]);
Rel_2_1_to_Rel_2_5_upgrade_30_date = datenum([2007,10,27]);

% ============ waiting for upgrade to Rel 2.1 =========================

disp(['The date of your scan is : ' num2str(scan_date) ' in [yyyy mm dd]'])

% convert the user given scan date to a days format
scan_date = datenum(scan_date);

% check the user given date against the known time points
if strcmpi(scanner, '1.5')
	if scan_date < VMS_to_windows_upgrade_15_date;
        date = '1.5_VMS_OS'; 
	else
        date = '1.5_Windows_OS'; 
    end
    
    % ============ waiting for upgrade to Rel 2.1 =========================
    %if scan_date < Rel_11_to_Rel_2_1_upgrade_15_date;
    release = 'Rel_11.x'; 
	%else
    %    release = 'Rel_2.1'; 
    %end
   % ============ waiting for upgrade to Rel 2.1 =========================
    
elseif strcmpi(scanner, '3.0')
    if scan_date < VMS_to_windows_upgrade_30_date
        date = '3.0_VMS_OS';
	else
        date = '3.0_Windows_OS';
    end
    
   % ============ waiting for upgrade to Rel 2.1 =========================
    if (scan_date < Rel_11_to_Rel_2_1_upgrade_30_date);
        release = 'Rel_11.x'; 
    elseif  ((Rel_11_to_Rel_2_1_upgrade_30_date < scan_date) & (scan_date < Rel_2_1_to_Rel_2_5_upgrade_30_date))
        release = 'Rel_2.1'; 
    elseif (Rel_2_1_to_Rel_2_5_upgrade_30_date < scan_date)
    	release = 'Rel_2.5';
    end
   % ============ waiting for upgrade to Rel 2.1 =========================
end

% Convert the angles in degrees into radians
ap = ap*2*pi/(360);
fh = fh*2*pi/(360);
rl = rl*2*pi/(360);

% Definition of Philips Master Gradient Table for YES-OVERPLUS and
% NO-Overplus, and definitions for jones30, invert-jones30, and jones30
% in the user defined option.  
 
if (strcmpi(release,'Rel_1.5') | strcmpi(release,'Rel_1.7') | strcmpi(release,'Rel_2.0') | strcmpi(release,'Rel_2.1') | strcmpi(release,'Rel_2.5'))  
	IDIFF_SCHEME_LOW_OVERPLUS =    [-0.4714,-0.9428,-0.9428;...
			     		-0.9428,-0.4714, 0.9428;...
			     		 0.9428,-0.9428, 0.4714;...
					-1.0000,-1.0000, 0.0000;...
					 0.0000,-1.0000, 1.0000;...
					 1.0000, 0.0000, 1.0000];


	IDIFF_SCHEME_MEDIUM_OVERPLUS = [-0.7071,-0.7071,-1.0000;...
					-0.7071,-0.7071, 1.0000;...
					 1.0000,-1.0000, 0.0000;...
					-0.1561,-0.9999,-0.9879;...
					 0.4091,-0.9894,-0.9240;...
					 0.8874,-0.4674,-0.9970;...
					 0.9297,-0.3866,-0.9930;...
					-0.9511,-0.7667,-0.7124;...
					 0.9954,-0.6945, 0.7259;...
					-0.9800,-0.3580, 0.9547;...
					-0.9992,-1.0000, 0.0392;...
					-0.3989,-0.9999, 0.9171;...
					 0.4082,-0.9923, 0.9213;...
					 0.9982,-0.9989, 0.0759;...
					 0.9919,-0.2899, 0.9655];


	if strcmpi(release,'Rel_2.5')
		% Philips changed one direction (11th from bottom) from [1.00000,-0.66820, 0.74400] to [-1.00000,-1.00000, 0.01110].
		
		IDIFF_SCHEME_HIGH_OVERPLUS =   [-0.70710,-0.70710,-1.00000;...
						-0.70710,-0.70710, 1.00000;...
						 1.00000,-1.00000, 0.00000;...
						-0.92390,-0.38270,-1.00000;...
						-0.29510,-0.95550,-1.00000;...
						 0.02780,-0.99960,-1.00000;...
						 0.59570,-0.80320,-1.00000;...
						 0.97570,-0.21910,-1.00000;...
						-0.92420,-0.38280,-0.99970;...
						-0.41420,-1.00000,-0.91020;...
						 0.41650,-0.99900,-0.91020;...
						 0.72830,-0.68740,-0.99850;...
						 1.00000,-0.41420,-0.91020;...
						-1.00000,-0.66820,-0.74400;...
						-0.66820,-1.00000,-0.74400;...
						 0.78560,-0.91070,-0.74400;...
						 1.00000,-0.66820,-0.74400;...
						-1.00000,-1.00000,-0.00030;...
						-1.00000,-0.66820, 0.74400;...
						 1.00000,-0.66820, 0.74400;...
						 0.66820,-1.00000, 0.74400;...
						-1.00000,-1.00000, 0.01110;...
						-0.90000,-0.60130, 0.91020;...
						-0.99850,-0.99850, 0.07740;...
						-0.41420,-1.00000, 0.91020;...
						 0.41420,-1.00000, 0.91020;...
						 1.00000,-1.00000, 0.01110;...
						 1.00000,-0.41420, 0.91020;...
						-0.99880,-0.99880, 0.06920;...
						 0.04910,-0.99880, 1.00000;...
						 0.99990,-0.99990, 0.01630;...
						 1.00000, 0.00000, 1.00000];
	else
		IDIFF_SCHEME_HIGH_OVERPLUS =   [-0.70710,-0.70710,-1.00000;...
						-0.70710,-0.70710, 1.00000;...
						 1.00000,-1.00000, 0.00000;...
						-0.92390,-0.38270,-1.00000;...
						-0.29510,-0.95550,-1.00000;...
						 0.02780,-0.99960,-1.00000;...
						 0.59570,-0.80320,-1.00000;...
						 0.97570,-0.21910,-1.00000;...
						-0.92420,-0.38280,-0.99970;...
						-0.41420,-1.00000,-0.91020;...
						 0.41650,-0.99900,-0.91020;...
						 0.72830,-0.68740,-0.99850;...
						 1.00000,-0.41420,-0.91020;...
						-1.00000,-0.66820,-0.74400;...
						-0.66820,-1.00000,-0.74400;...
						 0.78560,-0.91070,-0.74400;...
						 1.00000,-0.66820,-0.74400;...
						-1.00000,-1.00000,-0.00030;...
						-1.00000,-0.66820, 0.74400;...
						 1.00000,-0.66820, 0.74400;...
						 0.66820,-1.00000, 0.74400;...
						 1.00000,-0.66820, 0.74400;...
						-0.90000,-0.60130, 0.91020;...
						-0.99850,-0.99850, 0.07740;...
						-0.41420,-1.00000, 0.91020;...
						 0.41420,-1.00000, 0.91020;...
						 1.00000,-1.00000, 0.01110;...
						 1.00000,-0.41420, 0.91020;...
						-0.99880,-0.99880, 0.06920;...
						 0.04910,-0.99880, 1.00000;...
						 0.99990,-0.99990, 0.01630;...
						 1.00000, 0.00000, 1.00000];
	
	end

elseif ~(strcmpi(release,'Rel_1.5') | strcmpi(release,'Rel_1.7') | strcmpi(release,'Rel_2.0') | strcmpi(release,'Rel_2.1') |strcmpi(release,'Rel_2.5'))
	IDIFF_SCHEME_LOW_OVERPLUS = [0.9428,-0.4714,0.9428;...
                        	     0.4714,-0.9428,-0.9428;...
                        	     0.9428,0.9428,-0.4714;...
                        	     1.0,-1.0,0.0;...
                        	     1.0,0.0,-1.0;...
                        	     0.0,1.0,-1.0];


	IDIFF_SCHEME_MEDIUM_OVERPLUS = [0.7071,-0.7071,1.0000;...
                                	0.7071,-0.7071,-1.0000;...
                                	1.0000,1.0000,0.0000;...
                                	0.9999,-0.1561,0.9879;...
                                	0.9894,0.4091,0.9240;...
                                	0.4674,0.8874,0.9970;...
                                	0.3866,0.9297,0.9930;...
                                	0.7667,-0.9511,0.7124;...
                                	0.6945,0.9954,-0.7259;...
                                	0.3580,-0.9800,-0.9547;...
                                	1.0000,-0.9992,-0.0392;...
                                	0.9999,-0.3989,-0.9171;...
                                	0.9923,0.4082,-0.9213;...
                                	0.9989,0.9982,-0.0759;...
                                	0.2899,0.9919,-0.9655];

	IDIFF_SCHEME_HIGH_OVERPLUS = [0.3827,-0.9239,1.0000;...
                        	      0.9555,-0.2951,1.0000;...
                        	      0.9996,0.0278,1.0000;...
                        	      0.8032,0.5957,1.0000;...
                        	      0.2191,0.9757,1.0000;...
                        	      0.3828,-0.9242,0.9997;...
                        	      0.7071,-0.7071,1.0000;...
                        	      1.0000,-0.4142,0.9102;...
                        	      0.9990,0.4165,0.9102;...
                        	      0.6874,0.7283,0.9985;...
                        	      0.4142,1.0000,0.9102;...
                        	      0.6682,-1.0000,0.7440;...
                        	      1.0000,-0.6682,0.7440;...
                        	      0.9107,0.7856,0.7440;...
                        	      0.6682,1.0000,0.7440;...
                        	      1.0000,-1.0000,0.0003;...
                        	      1.0000,1.0000,0.0000;...
                        	      0.6682,-1.0000,-0.7440;...
                        	      0.6682,1.0000,-0.7440;...
                        	      1.0000,0.6682,-0.7440;...
                        	      0.6682,1.0000,-0.7440;...
                        	      0.6013,-0.9000,-0.9102;...
                        	      0.9985,-0.9985,-0.0774;...
                        	      1.0000,-0.4142,-0.9102;...
                        	      1.0000,0.4142,-0.9102;...
                        	      1.0000,1.0000,-0.0111;...
                        	      0.4142,1.0000,-0.9102;...
                        	      0.5624,-0.8269,-1.0000;...
                        	      0.9988,-0.9988,-0.0692;...
                        	      0.9988,0.0491,-1.0000;...
                        	      0.9999,0.9999,-0.0163;...
                        	      0.0000,1.0000,-1.0000];						
end	


IDIFF_SCHEME_LOW_NO_OVERPLUS = [1.0,0.0,0.0;...
                                0.0,1.0,0.0;...
                                0.0,0.0,1.0;...
                               -0.7044,-0.0881,-0.7044;...
                                0.7044,0.7044,0.0881;...
                                0.0881,0.7044,0.7044];
                        
                        
IDIFF_SCHEME_MEDIUM_NO_OVERPLUS = [1.0,0.0,0.0;...
                                   0.0,1.0,0.0;...
                                   0.0,0.0,1.0;...
                                  -0.1789,-0.1113,-0.9776;...
                                  -0.0635,0.3767,-0.9242;...
                                   0.710,0.0516,-0.7015;...
                                   0.6191,-0.4385,-0.6515;...  
                                   0.2424,0.7843,-0.5710;...
                                  -0.2589,-0.6180,-0.7423;...
                                  -0.8169,0.1697,-0.5513;...
                                  -0.8438,0.5261,-0.1060;...
                                  -0.2626,0.9548,-0.1389;...
                                   0.0001,0.9689,0.2476;...
                                   0.7453,0.6663,0.0242;...
                                   0.9726,0.2317,0.0209]; 

IDIFF_SCHEME_HIGH_NO_OVERPLUS = [1.0,0.0,0.0;...
                                 0.0,1.0,0.0;...
                                 0.0,0.0,1.0;...
                                -0.0424,-0.1146,-0.9925;...
                                 0.1749,-0.2005,-0.9639;...
                                 0.2323,-0.1626,-0.9590;...
                                 0.3675,0.0261,-0.9296;...
                                 0.1902,0.3744,-0.9076;...
                                -0.1168,0.8334,-0.5402;...
                                -0.2005,0.2527,-0.9466;...
                                -0.4958,0.1345,-0.8580;...
                                -0.0141,-0.6281,-0.7780;...
                                -0.7445,-0.1477,-0.6511;...
                                -0.7609,0.3204,-0.5643;...
                                -0.1809,0.9247,-0.3351;...
                                -0.6796,-0.4224,-0.5997;...
                                 0.7771,0.4707,-0.4178;...
                                 0.9242,-0.1036,-0.3677;...
                                 0.4685,-0.7674,-0.4378;...
                                 0.8817,-0.1893,-0.4322;...
                                 0.6904,0.7062,-0.1569;...
                                 0.2391,0.7571,-0.6080;...
                                -0.0578,0.9837,0.1703;...
                                -0.5368,0.8361,-0.1135;...
                                -0.9918,-0.1207,-0.0423;...
                                -0.9968,0.0709,-0.0379;...
                                -0.8724,0.4781,-0.1014;...
                                -0.2487,0.9335,0.2581;...
                                 0.1183,0.9919,-0.0471;...
                                 0.3376,0.8415,0.4218;...
                                 0.5286,0.8409,0.1163;...
                                 0.9969,0.0550,-0.0571]; 

IDIFF_SCHEME_JONES30 = [1,0,0;...
                        0.166,0.986,0;...             
                        -0.110,0.664,0.740;...   
                        0.901,-0.419,-0.110;...      
                        -0.169,-0.601, 0.781;...    
                        -0.815, -0.386, 0.433;...
                        0.656, 0.366, 0.660;...
                        0.582, 0.800, 0.143;...
                        0.900, 0.259, 0.350;...
                        0.693, -0.698, 0.178;...
                        0.357, -0.924, -0.140;...
                        0.543, -0.488, -0.683;...
                        -0.525, -0.396, 0.753;...
                        -0.639, 0.689, 0.341;...
                        -0.330, -0.013, -0.944;...
                        -0.524, -0.783, 0.335;...
                        0.609, -0.065, -0.791;...
                        0.220, -0.233, -0.947;...
                        -0.004, -0.910, -0.415;...
                        -0.511, 0.627, -0.589;...
                        0.414, 0.737, 0.535;...
                        -0.679, 0.139, -0.721;...
                        0.884, -0.296, 0.362;...
                        0.262, 0.432, 0.863;...
                        0.088, 0.185, -0.979;...
                        0.294, -0.907, 0.302;...
                        0.887, -0.089, -0.453;...
                        0.257, -0.443, 0.859;...
                        0.086, 0.867, -0.491;...
                        0.863, 0.504, -0.025];
                
IDIFF_SCHEME_INVERT_JONES30 = -IDIFF_SCHEME_JONES30;

IDIFF_SCHEME_USERDEFINED_JONES30 = IDIFF_SCHEME_JONES30;

% This line was the line actually used before the software upgrade.  NOTE, this is not a Jones30 direction,
% in fact, the first element has the WRONG sign by accident.  If data was
% collected with this direction, then we should use this direction for the
% analaysis. 
VMS_pre_upgrade_grad17 = [-0.609, -0.065, -0.791];
IDIFF_SCHEME_VMS_JONES30 = IDIFF_SCHEME_JONES30;
IDIFF_SCHEME_VMS_JONES30(17,:) = VMS_pre_upgrade_grad17;


% read in the supplied user defined or user-file gradient table.  It should have no more
% than 4 columns, with white space seperations.  ie.  0: 1.000 1.000 1.000 or 1.000 1.000 1.000
if (strcmpi(grad_choice,'user-defined') | strcmpi(grad_choice,'user-file'))
     disp(['This option is new, if it crashes, check the format of your supplied file']);
    IDIFF_SCHEME_USERDEFINED = load(supplied_grad_file);
    if size(IDIFF_SCHEME_USERDEFINED,2) > 4;
        error('The supplied gradient table file has more than 4 columns')
    elseif size(IDIFF_SCHEME_USERDEFINED,2) == 4;
        IDIFF_SCHEME_USERDEFINED = IDIFF_SCHEME_USERDEFINED(:,2:4);
    elseif size(IDIFF_SCHEME_USERDEFINED,2) == 3;

    end
    supplied_nvolumes = size(IDIFF_SCHEME_USERDEFINED,1);
end

if strcmpi(grad_choice,'user-defined')
    if ((nvolumes-1) == supplied_nvolumes)
        in = IDIFF_SCHEME_USERDEFINED;
        in_txt = 'IDIFF_SCHEME_USERDEFINED';
        disp('Based on the number of volumes in the par file, and the gradient table you provided')
        disp('The directions manually entered in the user defined field were used')
         space = 'LPH';
    else
        disp('You should not enter the 0,0,0 volume in the user-defined table you provide')
        error(['Your number of DWI volumes, ' num2str(nvolumes-1) ' in the data and the lines ' num2str(supplied_nvolumes) ' in your supplied user-defined gradient table file ' supplied_grad_file ' are in conflict'])
    end
elseif strcmpi(grad_choice,'user-file')
    if ((nvolumes-1) == supplied_nvolumes)
        in = IDIFF_SCHEME_USERDEFINED;
        in_txt = 'IDIFF_SCHEME_USERDEFINED';
        disp('Based on the number of volumes in the par file, and the gradient table you provided')
        disp('The directions entered via a text file and the user-file option were used')
         space = 'MPS';
    else
        disp('You should not enter the 0,0,0 volume in the user-file table you provide')
        error(['Your number of DWI volumes, ' num2str(nvolumes-1) ' in the data and the lines ' num2str(supplied_nvolumes) ' in your supplied user-file gradient table file' supplied_grad_file ' are in conflict'])
    end
elseif strcmpi(grad_choice,'jones30')
    if (nvolumes == 32 & (strcmpi(date,'1.5_Windows_OS') | strcmpi(date,'3.0_Windows_OS')))
        in = IDIFF_SCHEME_JONES30;
        in_txt = 'IDIFF_SCHEME_JONES30';
        disp('Based on the date of your scan and the number of volumes in the par file')
        disp('The Jones30 table from the windows OS was used');
        space = 'MPS';
    elseif (nvolumes == 35 & (strcmpi(date,'1.5_VMS_OS') | strcmpi(date,'3.0_VMS_OS')))
        in = IDIFF_SCHEME_VMS_JONES30;
        in_txt = 'IDIFF_SCHEME_VMS_JONES30';
        disp('Based on the date of your scan and the number of volumes in the par file')
        disp('The Jones30 table from the VMS OS was used');
        space = 'MPS';
    elseif (nvolumes == 31 & (strcmpi(date,'1.5_Windows_OS') | strcmpi(date,'3.0_Windows_OS')))
        in = IDIFF_SCHEME_USERDEFINED_JONES30;
        in_txt = 'IDIFF_SCHEME_USERDEFINED_JONES30'
        disp('Based on the date of your scan and the number of volumes in the par file')
        disp('The Jones30 table from the user defined field was used')
        space = 'LPH';
    else
        error(['Your number of volumes, ' num2str(nvolumes) ' scan date, and your input of grad_choice of ' grad_choice ' are in conflict'])
    end
elseif strcmpi(grad_choice,'invert-jones30')
    if (nvolumes == 32 & (strcmpi(date,'1.5_Windows_OS') | strcmpi(date,'3.0_Windows_OS')))
        in = IDIFF_SCHEME_INVERT_JONES30;
        in_txt = 'IDIFF_SCHEME_INVERT_JONES30';
        disp('Based on the date of your scan, your selection for grad_choice and the number of volumes in the par file')
        disp('The Invert-Jones30 table from the windows OS was used')
        space = 'MPS';
    else
        error(['Your number of volumes, ' num2str(nvolumes) ' scan date, and grad_choice of ' grad_choice ' are in conflict'])
    end 
elseif strcmpi(grad_choice,'yes-ovp-low')
    if nvolumes == 8 % assumes 6 directions + 1 b0 + 1 mean DWI
        in = IDIFF_SCHEME_LOW_OVERPLUS;
        in_txt = 'IDIFF_SCHEME_LOW_OVERPLUS';
        disp('Based on the information you provided and the number of volumes in the par file')
        disp('The low directional resolution, gradient overplus yes, table was used')
	if ~(strcmpi(release,'Rel_1.5') | strcmpi(release,'Rel_1.7') | strcmpi(release,'Rel_2.0') | strcmpi(release,'Rel_2.1') | strcmpi(release,'Rel_2.5'))
        	space = 'XYZ';
	elseif (strcmpi(release,'Rel_1.5') | strcmpi(release,'Rel_1.7') | strcmpi(release,'Rel_2.0') | strcmpi(release,'Rel_2.1') | strcmpi(release,'Rel_2.5'))
		space = 'LPH';
	end
	disp(['from release ' release ' which is defined in ' space ' space'])
    else
        error(['Your number of volumes, ' num2str(nvolumes) ' and grad_choice of ' grad_choice ' are in conflict'])       
    end
elseif strcmpi(grad_choice,'yes-ovp-medium')
    if nvolumes == 17 % assumes 15 directions + 1 b0 + 1 mean DWI
        in = IDIFF_SCHEME_MEDIUM_OVERPLUS;
        in_txt = 'IDIFF_SCHEME_MEDIUM_OVERPLUS';
	disp('Based on the information you provided and the number of volumes in the par file')
        disp('The medium directional resolution, gradient overplus yes, table was used')
        if ~(strcmpi(release,'Rel_1.5') | strcmpi(release,'Rel_1.7') | strcmpi(release,'Rel_2.0') | strcmpi(release,'Rel_2.1') | strcmpi(release,'Rel_2.5'))
        	space = 'XYZ';
	elseif (strcmpi(release,'Rel_1.5') | strcmpi(release,'Rel_1.7') | strcmpi(release,'Rel_2.0') | strcmpi(release,'Rel_2.1') | strcmpi(release,'Rel_2.5'))
		space = 'LPH';
	end
	disp(['from release ' release ' which is defined in ' space ' space'])
    else
       error(['Your number of volumes, ' num2str(nvolumes) ' and grad_choice of ' grad_choice ' are in conflict'])
    end
elseif strcmpi(grad_choice,'yes-ovp-high')
    if nvolumes == 34 % assumes 32 directions + 1 b0 + 1 mean DWI
        in = IDIFF_SCHEME_HIGH_OVERPLUS; 
        in_txt = 'IDIFF_SCHEME_HIGH_OVERPLUS';
        disp('Based on the information you provided and the number of volumes in the par file')
        disp('The high directional resolution, gradient overplus yes, table was used')
        if ~(strcmpi(release,'Rel_1.5') | strcmpi(release,'Rel_1.7') | strcmpi(release,'Rel_2.0') | strcmpi(release,'Rel_2.1') | strcmpi(release,'Rel_2.5'))
        	space = 'XYZ';
	elseif (strcmpi(release,'Rel_1.5') | strcmpi(release,'Rel_1.7') | strcmpi(release,'Rel_2.0') | strcmpi(release,'Rel_2.1') | strcmpi(release,'Rel_2.5'))
		space = 'LPH';
	end
	disp(['from release ' release ' which is defined in ' space ' space'])
    else
        error(['Your number of volumes, ' num2str(nvolumes) ' and grad_choice of ' grad_choice ' are in conflict'])
    end     
elseif strcmpi(grad_choice,'no-ovp-low')
    if nvolumes == 8
        in = IDIFF_SCHEME_LOW_NO_OVERPLUS;
        in_txt = 'IDIFF_SCHEME_LOW_NO_OVERPLUS';
        disp('Based on the date of your scan and the number of volumes in the par file')
        disp('The low directional resolution, gradient overplus no, table was used')
        space = 'MPS';
    else
        error(['Your number of volumes, ' num2str(nvolumes) ' and grad_choice of ' grad_choice ' are in conflict'])
    end
elseif strcmpi(grad_choice,'no-ovp-medium')
    if nvolumes == 17
        in = IDIFF_SCHEME_MEDIUM_NO_OVERPLUS;
        in_txt = 'IDIFF_SCHEME_MEDIUM_NO_OVERPLUS';
        disp('Based on the date of your scan and the number of volumes in the par file')
        disp('The medium directional resolution, gradient overplus no, table was used')
        space = 'MPS';
    else
        error(['Your number of volumes, ' num2str(nvolumes) ' and grad_choice of ' grad_choice ' are in conflict'])
    end
elseif strcmpi(grad_choice,'no-ovp-high')
    if nvolumes == 34
        in = IDIFF_SCHEME_HIGH_NO_OVERPLUS;
        in_txt = 'IDIFF_SCHEME_HIGH_NO_OVERPLUS';
        disp('Based on the date of your scan and the number of volumes in the par file')
        disp('The high directional resolution, gradient overplus no, table was used')
        space = 'MPS';
    else
        error(['Your number of volumes, ' num2str(nvolumes) ' and grad_choice of ' grad_choice ' are in conflict'])
    end
end

% Define some storage matrices
out = zeros(size(in));
rev_out = zeros(size(in));

% ==========================================================
% TRANSFORMATION DEFINITIONS 
% ==========================================================

% Transformations and reverse transformatins that we will use
% Definitions for these matrices were taken from Philips documentation

if strcmpi(patient_orientation,'sp')
    Tpo = [1,0,0;0,1,0;0,0,1];
    rev_Tpo = [1,0,0;0,1,0;0,0,1];
elseif strcmpi(patient_orientation,'pr')
    Tpo = [-1,0,0;0,-1,0;0,0,1];
    rev_Tpo = [-1,0,0;0,-1,0;0,0,1];  
elseif strcmpi(patient_orientation,'rd')
    Tpo = [0,-1,0;1,0,0;0,0,1];
    rev_Tpo = [0,1,0;-1,0,0;0,0,1];
elseif strcmpi(patient_orientation,'ld')
    Tpo = [0,1,0;-1,0,0;0,0,1];
    rev_Tpo = [0,-1,0;1,0,0;0,0,1];
end

if strcmpi(patient_position,'ff')
    Tpp = [0,-1,0;-1,0,0;0,0,1];
    rev_Tpp = [0,-1,0;-1,0,0;0,0,-1];
elseif strcmpi(patient_position,'hf')
    Tpp = [0,1,0;-1,0,0;0,0,-1];
    rev_Tpp = [0,-1,0;1,0,0;0,0,-1];
end

Tpom = Tpo*Tpp;
rev_Tpom = rev_Tpp*rev_Tpo;
    
% Definitions for Trl,Tap,Tfh and Tang
Trl = [1,0,0; 0, cos(rl), -sin(rl);0,sin(rl),cos(rl)];
Tap = [cos(ap),0,sin(ap); 0,1,0; -sin(ap),0,cos(ap)];
Tfh = [cos(fh),-sin(fh),0; sin(fh),cos(fh),0; 0,0,1];
Tang = Trl*Tap*Tfh;

rev_Trl = [1,0,0; 0, cos(rl), sin(rl);0,-sin(rl),cos(rl)];
rev_Tap = [cos(ap),0,-sin(ap); 0,1,0; sin(ap),0,cos(ap)];
rev_Tfh = [cos(fh),sin(fh),0;-sin(fh),cos(fh),0;0,0,1];
rev_Tang = rev_Tfh*rev_Tap*rev_Trl;

% Definitions for Tsom
if strcmpi(slice_orientation,'sag')
    Tsom = [0,0,-1;0,-1,0;1,0,0];
    rev_Tsom = [0,0,1;0,-1,0;-1,0,0];
elseif strcmpi(slice_orientation,'cor')
    Tsom = [0,-1,0;0,0,1;1,0,0];
    rev_Tsom = [0,0,1;-1,0,0;0,1,0];
elseif strcmpi(slice_orientation,'tra')
    Tsom = [0,-1,0; -1,0,0; 0,0,1];
    rev_Tsom = [0,-1,0;-1,0,0;0,0,1];
end

% Definitions for Tprep_par Tprep_per & Tfsd_m, Tfsd_p, Tfsd_s

Tprep_par = [1,0,0;0,1,0;0,0,1];
rev_Tprep_par = [1,0,0;0,1,0;0,0,1];
Tprep_per = [0,-1,0;1,0,0;0,0,1];
rev_Tprep_per = [0,1,0;-1,0,0;0,0,1];

Tfsd_m = [-1,0,0;0,1,0;0,0,1];
rev_Tfsd_m = [-1,0,0;0,1,0;0,0,1];
Tfsd_p = [1,0,0;0,-1,0;0,0,1];
rev_Tfsd_p = [1,0,0;0,-1,0;0,0,1];
Tfsd_s = [1,0,0;0,1,0;0,0,-1];
rev_Tfsd_s = [1,0,0;0,1,0;0,0,-1];

if strcmpi(slice_orientation,'tra')
    if strcmpi(foldover,'AP')
        Tprep = Tprep_per;
        rev_Tprep = rev_Tprep_per;
        if strcmpi(fat_shift,'A')
            Tfsd = Tfsd_m;
            rev_Tfsd = rev_Tfsd_m;
        elseif strcmpi(fat_shift,'P')
            Tfsd = Tfsd_p;
            rev_Tfsd = rev_Tfsd_p;
        end
    elseif strcmpi(foldover,'RL')
        Tprep = Tprep_par;
        rev_Tprep = rev_Tprep_par;
        if strcmpi(fat_shift,'R')
            Tfsd = Tfsd_p;
            rev_Tfsd = rev_Tfsd_p;
        elseif strcmpi(fat_shift,'L')
            Tfsd = Tfsd_m;
            rev_Tfsd = rev_Tfsd_m;
        end
    end
    
elseif strcmpi(slice_orientation,'cor')
    if strcmpi(foldover,'FH')
        Tprep = Tprep_per;
        rev_Tprep = rev_Tprep_per;
        if strcmpi(fat_shift,'F')
            Tfsd = Tfsd_p;
            rev_Tfsd = rev_Tfsd_p;
        elseif strcmpi(fat_shift,'H')
            Tfsd = Tfsd_m;
            rev_Tfsd = rev_Tfsd_m;
        end
    elseif strcmpi(foldover,'RL')
        Tprep = Tprep_par;
        rev_Tprep = rev_Tprep_par;
        if strcmpi(fat_shift,'R')
            Tfsd = Tfsd_p;
            rev_Tfsd = rev_Tfsd_p;
        elseif strcmpi(fat_shift,'L')
            Tfsd = Tfsd_m;
            rev_Tfsd = rev_Tfsd_m;
        end
    end
    
elseif strcmpi(slice_orientation,'sag')
    if strcmpi(foldover,'FH')
        Tprep = Tprep_per;
        rev_Tprep = rev_Tprep_per;
        if strcmpi(fat_shift,'F')
            Tfsd = Tfsd_p;
            rev_Tfsd = rev_Tfsd_p;
        elseif strcmpi(fat_shift,'H')
            Tfsd = Tfsd_m;
            rev_Tfsd = rev_Tfsd_m;
        end
    elseif strcmpi(foldover,'AP')
        Tprep = Tprep_par;
        rev_Tprep = rev_Tprep_par;
        if strcmpi(fat_shift,'A')
            Tfsd = Tfsd_p;
            rev_Tfsd = rev_Tfsd_p;
        elseif strcmpi(fat_shift,'P')
            Tfsd = Tfsd_m;
            rev_Tfsd = rev_Tfsd_m;
        end
    end
end

% ==========================================
% END OF PHILIPS TRANSFORMATION DEFINITIONS
% ==========================================

% Transformation needed to go from Philips NWV coordinate space to DTIstudio
% coordinate space
DTIextra = [0,-1,0;-1,0,0;0,0,1];
rev_DTIextra = [0,-1,0;-1,0,0;0,0,1];

% ======================================
% APPLICATION OF THE TRANSFORMATIONS
% ======================================
if strcmpi(space,'LPH')
    for i = 1:length(in)
        if strcmpi(DTI_studio,'n')
            out(i,:) = (rev_Tsom*rev_Tang*in(i,:)')';
            message = 'Ang_corrected_table is a NVW pixel space compatible table from LPH origin';
            % to check things, apply the reverse transformation to go
            % from out back to in
            rev_out(i,:) = (Tang*Tsom*out(i,:)')';
        elseif strcmpi(DTI_studio,'y')
            out(i,:) = (DTIextra*rev_Tsom*rev_Tang*in(i,:)')';
            message = 'Ang_corrected_table is a DTIstudio compatible table from LPH origin';
            % to check things, apply the reverse transformation to go
            % from out back to in
            rev_out(i,:) = (Tang*Tsom*rev_DTIextra*out(i,:)')';
        end
	end
    %disp(message)
elseif strcmpi(space,'XYZ')
	for i = 1:length(in)
        if strcmpi(DTI_studio,'n')
            out(i,:) = (rev_Tsom*rev_Tang*Tpom*in(i,:)')';
            message = 'Ang_corrected_table is a NVW pixel space compatible table from XYZ origin';
            % to check things, apply the reverse transformation to go
            % from out back to in
            rev_out(i,:) = (rev_Tpom*Tang*Tsom*out(i,:)')';
        elseif strcmpi(DTI_studio,'y')
            out(i,:) = (DTIextra*rev_Tsom*rev_Tang*Tpom*in(i,:)')';
            message = 'Ang_corrected_table is a DTIstudio compatible table from XYZ origin ';
            % to check things, apply the reverse transformation to go
            % from out back to in
            rev_out(i,:) = (rev_Tpom*Tang*Tsom*rev_DTIextra*out(i,:)')';
        end
    end
    %disp(message)
elseif strcmpi(space,'MPS')
    for i = 1:length(in)
        if strcmpi(DTI_studio,'n')
            out(i,:) = (Tprep*Tfsd*in(i,:)')';
            message = 'Ang_corrected_table is a NVW pixel space compatible table from MPS origin';
            % to check things, apply the reverse transformation to go
            % from out back to in
            rev_out(i,:) = (rev_Tfsd*rev_Tprep*out(i,:)')';
        elseif strcmpi(DTI_studio,'y')
            out(i,:) = (DTIextra*Tprep*Tfsd*in(i,:)')';
            message = 'Ang_corrected_table is a DTIstudio compatible table from MPS origin';
            % to check things, apply the reverse transformation to go
            % from out back to in
            rev_out(i,:) = (rev_Tfsd*rev_Tprep*rev_DTIextra*out(i,:)')';
        end
    end
    %disp(message)
end

% Normalize the non zero vectors
for ii = 1:length(out)
	if ((out(ii,1) == 0) & (out(ii,2) == 0) & (out(ii,3) == 0))
		out(ii,:) = out(ii,:);
	elseif ((out(ii,1) == 100) & (out(ii,2) == 100) & (out(ii,3) == 100))
		out(ii,:) = out(ii,:);
	else
		out(ii,:) = (1/norm(out(ii,:)))*out(ii,:);
	end
end

% Before the upgrade on the 1.5T, 5 b=0 images and no mean DWI were the output
if strcmpi(in_txt,'IDIFF_SCHEME_VMS_JONES30')
    if strcmpi(sort_images,'n')
        error('Sorry, I have no idea how VMS worked with sort images = no')
    elseif strcmpi(sort_images,'y')
        out = [0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;out];
        % the b=0 volume is first complete volume, 
    end
% The User Defined option does not produce a mean DWI volume in the rec file    
elseif strcmpi(in_txt,'IDIFF_SCHEME_USERDEFINED_JONES30')
    if strcmpi(sort_images,'n')
        out = [0,0,0;out];
        % the b=0 volume is the first sub volume,
    elseif strcmpi(sort_images,'y')
        out = [0,0,0;out];
        % the b=0 volume is first complete volume, 
    end 
% The User Defined option does not produce a mean DWI volume in the rec file    
elseif strcmpi(in_txt,'IDIFF_SCHEME_USERDEFINED')
    if strcmpi(sort_images,'n')
        if (strcmpi(par_version,'4.1') | strcmpi(par_version,'4.2'))
            out = [0,0,0; out]; 
            % the b=0 volume is first sub volume and must be listed first
            % to work with DTIstudio
        elseif strcmpi(par_version,'4')
            out = [0,0,0; out]; 
            % the b=0 volume is first sub volume and must be listed first
            % to work with DTIstudio
        elseif strcmpi(par_version,'3')
            error('Sorry, I have no idea how V3 worked with sort images = no')
        else
            error('unsupported and unknown par file vesion, must be 3, 4, 4.1 or 4.2')
        end      
    elseif strcmpi(sort_images,'y')
        if (strcmpi(par_version,'4.1') | strcmpi(par_version,'4.2'))
            out = [out; 0,0,0]; 
            % the b=0 volume is the last volume
        elseif (strcmpi(par_version,'3') | strcmpi(par_version,'4'))
            out = [0,0,0; out]; 
            % the b=0 volume is the first volume
        else
           error('unsupported and unknown par file vesion, must be 3, 4, 4.1 or 4.2')
        end     
    end
% All The Philips Tables and the user-file produce a mean DWI volume in the rec file 
else
    if strcmpi(sort_images,'n')
        if (strcmpi(par_version,'4.1') | strcmpi(par_version,'4.2'))
            out = [0,0,0; out; 100,100,100]; 
            % the b=0 volume is first sub volume and must be listed first
            % to work with DTIstudio
        elseif strcmpi(par_version,'4')
             out = [0,0,0; out; 100,100,100]; 
            % the b=0 volume is first sub volume and must be listed first
            % to work with DTIstudio
        elseif strcmpi(par_version,'3')
            error('Sorry, I have no idea how V3 worked with sort images = no')
        else
            error('unsupported and unknown par file vesion, must be 3, 4, 4.1 or 4.2')
        end
    elseif strcmpi(sort_images,'y')
         if (strcmpi(par_version,'4.1') | strcmpi(par_version,'4.2'))
            out = [out; 0,0,0; 100,100,100]; 
            % the b=0 volume is the 2nd last complete volume
        elseif (strcmpi(par_version,'3') | strcmpi(par_version,'4'))
            out = [0,0,0; out; 100,100,100]; 
            % the b=0 volume is the first complete volume
        else
            error('unsupported and unknown par file vesion, must be 3, 4, 4.1 or 4.2')
        end
    end
end

input_table = in;
Ang_corrected_table = out;

% =========================================================
% ========================================================
% END OF ANGULATION CORRECTION
% ========================================================
% ========================================================


% =========================================================
% ========================================================
% START OF REGISTRATION CORRECTION
% ========================================================
% ========================================================
function [T_and_Ang_corrected_table] = registration_correction(xfm_air_directory, name, Ang_corrected_table);

% August 2, 2006 | added fullfile and strtrim
% August 30, 2006 | added check for slash on xfm_air_direcotry


% add check to make sure xfm_air_directory has a slash at the end
if strcmpi(xfm_air_directory(end),'/')

else
	xfm_air_directory = [xfm_air_directory '/'];
end


disp(['Looking in ' xfm_air_directory ' for ' name '*.xfm or ' name '*.air files'])
% If there are .xfm files in the folder
if ~isempty(dir([xfm_air_directory name '*.xfm']))
    B = dir([xfm_air_directory name '*.xfm']);
	xfm_files = char(B.name);
	xfm_files = unique(xfm_files,'rows');
	
    % check that the gradient dynamics are ordered in 1 to n order
    % so that they can be used to transform the matching line in the
    % gradient table.  Note that issorted checks ASCII order. As the
    % xfm_files should only differ in the dynamic number, this should work.
    
    if ~issorted(xfm_files,'rows')
        disp('files were not sorted...so sorting was done')
        xfm_files = sortrows(xfm_files);
    end
    
	transforms = zeros([size(xfm_files,1) 16]);
    
	for ii = 1:size(xfm_files,1)
        	%[xfm_files(ii,:)]
		transforms(ii,:) = textread(fullfile(xfm_air_directory, strtrim(xfm_files(ii,:))),'%f')';
	end
    
    xfm_air_files = xfm_files;
    
elseif ~isempty(dir([xfm_air_directory name '*.air']))
    B = dir([xfm_air_directory name '*.air']);
	air_files = char(B.name);
	air_files = unique(air_files,'rows');
    
    % check that the gradient dynamics are ordered in 1 to n order
    % so that they can be used to transform the matching line in the
    % gradient table.  Note that issorted checks ASCII order. As the
    % xfm_files should only differ in the dynamic number, this should work.
    
    if ~issorted(air_files,'rows')
        disp('files were not sorted...so sorting was done')
        air_files = sortrows(air_files);
    end
    
    transforms = zeros([size(air_files,1) 16]);
    
    for ii = 1:size(air_files,1)
        [air_files(ii,:)]
	matrix = scanair_internal(fullfile(xfm_air_directory, strtrim(air_files(ii,:))));
        transforms(ii,:) = matrix(:)';
    end 
    xfm_air_files = air_files;
end

% Check the dimensions of the Ang_corrected_table and transforms
if size(Ang_corrected_table,1) ~= size(transforms,1)
    disp(['The size of the gradient table is' num2str(size(Ang_corrected_table,1)) 'by' num2str(size(Ang_corrected_table,2))])
    disp(['The size of the gradient table is' num2str(size(transforms,1)) 'by' num2str(size(transforms,2))])
    error('The number of rows in the transform matrix and the gradient table do not match')     
end

% need a line to make sure we don't normalize the mean DWI volume which we
% want to label with 100,100,100.
Ang_corrected_table
for j = 1:size(transforms,1)
    % reshape to a 4 by 4 matrix and take the Transpose
    T = reshape(transforms(j,:),[4 4])';
    % Only keep the rotations, disgard the translations
    T = T(1:3,1:3);
    if (all(Ang_corrected_table(j,:) == 100))
        T_and_Ang_corrected_table(j,:) = [100,100,100];
    end
    
    if all(~Ang_corrected_table(j,:))
        T_and_Ang_corrected_table(j,:) = [0,0,0];
    end
    
    if (~(all(Ang_corrected_table(j,:) == 100)) & (~all(~Ang_corrected_table(j,:))))
    
        T_and_Ang_corrected_table(j,:) = (T*Ang_corrected_table(j,:)')';
        % normalize the new direction
        if ~((T_and_Ang_corrected_table(j,1) == 0) & (T_and_Ang_corrected_table(j,2) == 0) & (T_and_Ang_corrected_table(j,3)==0))
            T_and_Ang_corrected_table(j,:) = (1/norm(T_and_Ang_corrected_table(j,:)))*T_and_Ang_corrected_table(j,:);
        end
    end
end
    
% =========================================================
% ========================================================
% END OF REGISTRATION CORRECTION
% ========================================================
% ========================================================

function [nrows, ncols, nslices, nechoes, ndynamics, nphases, A, header, par_version] = int_getPARinfo(filename_par)
%  addapted from Craig Jones's parseHeader - Parse a Philips PAR file for some important parameters
%

%  Craig Jones (craig@mri.jhu.edu)  
%  20040616 - fixed up reading in the header for stupid Windoze boxes
%  Jonathan Farrell (JonFarrell@jhu.edu)
%  20050628 - added the examination date, patient_position, and preparation_direction to the output header
%  20061213 - added code to read in V4.1 par files
%  20070411 - fixed code to return 3, 4 or 4.1 for par version
%  20071113 - fixed parsing statements for #  sl in v4.2 par files

nrows = 0; ncols = 0; nslices = 0; nechoes = 0; ndynamics = 0; A = []; header = []; par_version = 0;

header.filename = filename_par;

line = '';
fp = fopen(filename_par, 'rt');

if( fp == -1 )
	if( isunix )
		person = getenv('USER');
	elseif( ispc )
		person = getenv('USERNAME');
	else
		person = 'matlab user';
	end

	disp(sprintf('readrec:  I''m sorry, %s, the file %s does not exist', getenv('USER'), filename_par));
	return;
end

firstheader = 1;
par_version = 0;

%line = fgetl(fp);
while( 1 )
    line = fgetl(fp);

    if ((strncmp('#sl', line, 3) == 1) | (strncmp('# sl', line, 3) == 1)| (strncmp('#  sl ', line, 6))), break, end;
    
    if( firstheader & line(1) == '#' ) 
        % look for V4.1 or something similar
		[s,f] = regexp(line, '\w*[vV][0-9]\.\w*');
		if( ~isempty(s) )
            par_version = str2num( line(s(end)+1:f(end)) );
        else
            % look for V4 or V3 or something similar
            [s,f] = regexp(line, '\w*[vV][0-9]\w*');
            if( ~isempty(s) )
                par_version = str2num( line(s(end)+1:f(end)) );
            end
        end
	end
    
	%% We are not the in first header any more
    if ( strncmp('#', line, 1) ~= 1) 
		firstheader = 0;
	end
    
    %%  Look for off center values
    if( findstr('Off Centre ', line) > 0 )
        aa = line(findstr(':', line)+1:end);
       	header.off_center = str2num(aa);
		header.annotation = 'ap fh rl';
    end

    %%  Look for off center angulation
    if( findstr('Angulation ', line) > 0 )
        aa = line(findstr(':', line)+1:end);
       	header.angulation = str2num(aa);
    end

    %%  Look for the FOV
    if( findstr('FOV (', line) > 0 )
        aa = line(findstr(':', line)+1:end);
       	header.fov = str2num(aa);
    end

    %%  Look for number of rows and columns
    if( findstr('Recon resolution', line) > 0 )
        aa = line(findstr(':', line)+1:end);
        aa = str2num(aa);
        nrows = aa(1);  ncols = aa(2);
		header.recon_resolution = aa;
    end

    %%  Look for number of slices
    if( findstr('number of slices', line) > 0 )
        aa = line(findstr(':', line)+1:end);
        aa = str2num(aa);
        nslices = aa(1);
    end

    %%  Look for number of slices
    if( findstr('number of cardiac phases', line) > 0 )
        aa = line(findstr(':', line)+1:end);
        aa = str2num(aa);
        nphases = aa(1);
    end
    
    % ---------------------
    % added by Jonathan Farrell
    %%  Look for the Examination date
    if( findstr('Examination date/time', line) > 0 )
        aa = line(findstr(':', line)+1:end);
        aa = strrep(aa,'.',' ');
        aa = strrep(aa,'/',' ');
        aa = str2num(aa);
        header.examination_date = [aa(1) aa(2) aa(3)];
        header.date_annotation = '[yyyy mm dd]';
    end
    
    if par_version > 3
        if ( findstr('Patient position', line) > 0 )
            aa = line(findstr(':', line)+1:end);
            header.patient_position = aa;
        end
        
        if ( findstr('Preparation direction', line) > 0 )
            aa = line(findstr(':', line)+1:end);
            header.preparation_direction = aa;
        end  
    end
    
    % ---------------------------
end

line = fgetl(fp);

A = [];
ii = 1;
while( 1 )

    line = fgetl(fp);
    if( length(line) < 2 ), break; end
    
    A(ii,:) = str2num(line);
    ii = ii + 1;
end

ndynamics = length(unique(A(:,3)));
nechoes = length(unique(A(:,2)));

%  Added for the new PAR files.
if( nrows == 0 )
	nrows = A(1,10);
	ncols = A(1,11);
end

fclose(fp);

% =========================================================
% ========================================================
% START OF .AIR FILE READ IN FUNCTION
% ========================================================
% ========================================================

function [matrix] = scanair_internal(airfile)

% ==========================================
% HEADER
% ==========================================

% Name: scanair_internal.m

% Author:   Jonathan Farrell ( JonFarrell@jhu.edu )
%           F.M. Kirby Research Center for Functional Brain Imaging
%           Room G25
%           Kennedy Krieger Institute
%           Baltimore, MD 21218

% Creation Date:  October 10, 2005

% History of Updates:


% ========================================
% PURPOSE:  
% ========================================

% to mimic the scanair function provided by AIR, but do it in matlab to
% make it OS system (unix vs windows) independent

% The equivilant code using the AIR scanair package would be
% airfile = '/home/jfarrell/matlab_jon/img_1_0.air'
% [d,A] = unix(sprintf(['/usr/local/air5.2.5/scanair_16u ' airfile ' -r']));
% matrix = str2num(A(strfind(A,'[')+1:strfind(A,']')-1))       
% transforms(1,:) = matrix(:)'

% ========================================
% BEGIN CODE:  
% ========================================
endian = 'l';
fid = fopen(airfile,'rb',endian);
e = fread(fid,[4 4],'double');
s_file = fgets(fid, 128);
s_bits = fread(fid,1,'uint32'); % This is the number of bits per pixel in the dataset

% if the bits are not reasonable, try the other endian format
if (s_bits ~= 16)
    disp('Bits per pixel for little endian is not 16, ... Trying other endian')
    fclose(fid);
    endian = setdiff(['l','b'],endian); 
    fid = fopen(airfile,'rb',endian);
    e = fread(fid,[4 4],'double');
    s_file = fgets(fid, 128);
    s_bits = fread(fid,1,'uint32'); % This is the number of bits per pixel in the dataset
    if (s_bits ~= 16)
        error('Tried little and big endian, could not match bits per pixel of 16 critea')
    else
        disp(['Used ' endian ' endian format'])
    end
end

s_xdim = fread(fid,1,'uint32');
s_ydim = fread(fid,1,'uint32');
s_zdim = fread(fid,1,'uint32');
s_xsize = fread(fid,1,'double');
s_ysize = fread(fid,1,'double');
s_zsize = fread(fid,1,'double'); 

r_file = fgets(fid, 128);
r_bits = fread(fid,1,'uint32'); % This is the number of bits per pixel in the dataset
r_xdim = fread(fid,1,'uint32');
r_ydim = fread(fid,1,'uint32');
r_zdim = fread(fid,1,'uint32');
r_xsize = fread(fid,1,'double');
r_ysize = fread(fid,1,'double');
r_zsize = fread(fid,1,'double'); 

% to put the .air matrix in real world units in mm
% this should give the same result as the 'r' option when using the AIR
% based scanair function.  The real world transformation coordinates in mm
% make the most sense when thinking about how to correct the gradient
% tables.

pixel_size_s = min([s_xsize,s_ysize,s_zsize]);
matrix = e;
for j = 1:3
    for i = 1:4
        matrix(i,j) = matrix(i,j)/pixel_size_s;
    end
end
for j = 1:4
    matrix(1,j) = matrix(1,j)*r_xsize;
    matrix(2,j) = matrix(2,j)*r_ysize;
    matrix(3,j) = matrix(3,j)*r_zsize;
end

% =========================================================
% ========================================================
% END OF .AIR FILE READ IN FUNCTION
% ========================================================
% ========================================================
