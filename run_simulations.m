% --- Setup experiment params ----- %
phase_lag = pi / 2;
n_trials = 100;
InducedScale = 0.5; % a.k.a. signal to noise ratio
% ----
% don't change this section
gain_svd_th = 0.001;
EvokedScale = 0;
lambda = 100;
is_imag = false;
n_connections = 300;
detection_diam = 0.015;
n_steps = 100;
NPI = [1,2,3];
% ---------------------------------- %



[HM, CT, Trials, Ctx, XYZGenAct] = SimulateData(phase_lag, n_trials, gain_svd_th,...
                                                InducedScale, EvokedScale);


CT_reshape = reshape(mean(CT, 2), sqrt(size(CT,1)), sqrt(size(CT,1)));


[A, Ps, Cs, IND] = DICS(CT_reshape, HM.gain, lambda, is_imag);

con_inds = threshold_connections(Cs, n_connections, IND);
drawConnectionsOnBrain(con_inds, HM.GridLoc, Ctx);


 [SPC, TPR, PPV] = GenerateScores(Cs, detection_diam, HM.GridLoc, IND,...
                               n_steps, XYZGenAct, NPI);

 % ---- Plot precision-recall curve ---- %
 figure;
 plot(TPR, PPV);
