clear all;
close all;
	
%-------------------------------------------------------
% 1. Baseline Calibration
%-------------------------------------------------------

gamma = 0;

save paramfile gamma;

% sets we loop over
gamma_set = [0,.1,.2,.3,.4,.5,.6,.7,.8,.9];

%-------------------------------------------------------
% 2. Looping over Dynare code: rule-of-thumb fraction
%-------------------------------------------------------

for gamma = gamma_set

	save paramfile gamma;

    dynare nk_rebates_mpc_monthly_loops;

%   The sheet name must be a positive integer, so I had to transform gamma

    sheetgamma = gamma*10 + 1;
    
    % xlswrite does not work on mac / linux, replaced with writematrix
%     xlswrite('nk_rebate_loop_irfs.xls', A, sheetgamma)
    writematrix(A,'baseline_irfs.xlsx','Sheet',sheetgamma)

% Either suspend dropbox sync or set pause time to avoid the "terminate called after throwing ..." error
%    pause(4)
    
end



