function options = gpmp_optimset_init()
%GPMP_OPTIMSET_INIT Get a default struct of gpmp optimset options

options.MaxIter = 100;
options.MinIter = 1;
options.TolFun = 1e-3;
options.TolCon = 1e-3;
options.RecCostHis = true; %record cost history
options.RecStateHis = true; %record variable states history
options.FixedStateId = []; %fixed support states indexes
options.UseOSQP = false; %use osqp solver flag
end