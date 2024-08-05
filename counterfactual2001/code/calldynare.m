%% This .m file calls dynare and exports the results

dynare tank_rebates_2001_consadj.mod

A = [oo_.exo_simul, y, c, co, cr, t, to, tr, h, ho, hr, w, iv, u, rk, r, pi, k, b] ;

csvwrite('../output/baseline_2001_irfs.csv',A);

%Delete all files dynare just created
delete tank_rebates_2001_consadj_dynamic.m
delete tank_rebates_2001_consadj_results.mat
delete tank_rebates_2001_consadj_set_auxiliary_variables.m
delete tank_rebates_2001_consadj_static.m
delete tank_rebates_2001_consadj_steadystate2.m
delete tank_rebates_2001_consadj.log
delete tank_rebates_2001_consadj.m
rmdir('tank_rebates_2001_consadj','s')


exit
