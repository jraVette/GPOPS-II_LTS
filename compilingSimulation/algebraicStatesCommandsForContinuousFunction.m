%This function is appened to the continuous function in
%calculateAlgebraicStates.m to return the values of the algebraic states.
%If you include them natively, it slows the optimizaiton time
%significantly by the size the of the algebraic states saved.

phaseout.algebraicStates.sa_L1.meas = alpha_L1;
phaseout.algebraicStates.sa_R1.meas = alpha_R1;
phaseout.algebraicStates.sa_L2.meas = alpha_L2;
phaseout.algebraicStates.sa_R2.meas = alpha_R2;

phaseout.algebraicStates.slipRatio_L1.meas = kappa_L1;
phaseout.algebraicStates.slipRatio_R1.meas = kappa_R1;
phaseout.algebraicStates.slipRatio_L2.meas = kappa_L2;
phaseout.algebraicStates.slipRatio_R2.meas = kappa_R2;

phaseout.algebraicStates.fx_L1.meas = fx_L1;
phaseout.algebraicStates.fx_R1.meas = fx_R1;
phaseout.algebraicStates.fx_L2.meas = fx_L2;
phaseout.algebraicStates.fx_R2.meas = fx_R2;

phaseout.algebraicStates.fy_L1.meas = fy_L1;
phaseout.algebraicStates.fy_R1.meas = fy_R1;
phaseout.algebraicStates.fy_L2.meas = fy_L2;
phaseout.algebraicStates.fy_R2.meas = fy_R2;

phaseout.algebraicStates.muX_L1.meas = muX_L1;
phaseout.algebraicStates.muX_R1.meas = muX_R1;
phaseout.algebraicStates.muX_L2.meas = muX_L2;
phaseout.algebraicStates.muX_R2.meas = muX_R2;

phaseout.algebraicStates.muY_L1.meas = muY_L1;
phaseout.algebraicStates.muY_R1.meas = muY_R1;
phaseout.algebraicStates.muY_L2.meas = muY_L2;
phaseout.algebraicStates.muY_R2.meas = muY_R2;

% phaseout.algebraicStates.FxMax_L1.meas = FxMax_L1;
% phaseout.algebraicStates.FxMax_R1.meas = FxMax_R1;
% phaseout.algebraicStates.FxMax_L2.meas = FxMax_L2;
% phaseout.algebraicStates.FxMax_R2.meas = FxMax_R2;
% 
% phaseout.algebraicStates.FyMax_L1.meas = FyMax_L1;
% phaseout.algebraicStates.FyMax_R1.meas = FyMax_R1;
% phaseout.algebraicStates.FyMax_L2.meas = FyMax_L2;
% phaseout.algebraicStates.FyMax_R2.meas = FyMax_R2;

phaseout.algebraicStates.kappa_n_L1.meas = kappa_n_L1;
phaseout.algebraicStates.kappa_n_R1.meas = kappa_n_R1;
phaseout.algebraicStates.kappa_n_L2.meas = kappa_n_L2;
phaseout.algebraicStates.kappa_n_R2.meas = kappa_n_R2;

% phaseout.algebraicStates.alpha_n_L1.meas = alpha_n_L1;
% phaseout.algebraicStates.alpha_n_R1.meas = alpha_n_R1;
% phaseout.algebraicStates.alpha_n_L2.meas = alpha_n_L2;
% phaseout.algebraicStates.alpha_n_R2.meas = alpha_n_R2;

phaseout.algebraicStates.rho_L1.meas = rho_L1;
phaseout.algebraicStates.rho_R1.meas = rho_R1;
phaseout.algebraicStates.rho_L2.meas = rho_L2;
phaseout.algebraicStates.rho_R2.meas = rho_R2;

phaseout.algebraicStates.eff_L1.meas = eff_L1;
phaseout.algebraicStates.eff_R1.meas = eff_R1;
phaseout.algebraicStates.eff_L2.meas = eff_L2;
phaseout.algebraicStates.eff_R2.meas = eff_R2;


phaseout.algebraicStates.FX.meas = FX;
phaseout.algebraicStates.FY.meas = FY;
phaseout.algebraicStates.ax.meas = FX./m;
phaseout.algebraicStates.ay.meas = FY./m;

phaseout.algebraicStates.T_drive_L1.meas = T_drive_L1;
phaseout.algebraicStates.T_drive_R1.meas = T_drive_R1;
phaseout.algebraicStates.T_drive_L2.meas = T_drive_L2;
phaseout.algebraicStates.T_drive_R2.meas = T_drive_R2;

phaseout.algebraicStates.kt.meas = kt;
phaseout.algebraicStates.diffTorqueTransfer.meas = kd*(omega_L2 - omega_R2);
phaseout.algebraicStates.tPlus.meas =tPlus;
phaseout.algebraicStates.tMinus.meas =(1-tPlus);

phaseout.algebraicStates.enginePercent.meas = percentEnginePowerUsed;

% % phaseout.algebraicStates.vx_scaled.meas = vx;
% % phaseout.algebraicStates.vy_scaled.meas = vy;
% % phaseout.algebraicStates.r_scaled.meas = r;
% % phaseout.algebraicStates.omega_L1_scaled.meas = omega_L1;
% % phaseout.algebraicStates.omega_R1_scaled.meas = omega_R1;
% % phaseout.algebraicStates.omega_L2_scaled.meas = omega_L2;
% % phaseout.algebraicStates.omega_R2_scaled.meas = omega_R2;
% % phaseout.algebraicStates.T_scaled.meas = T;
% % phaseout.algebraicStates.ey_scaled.meas = ey;
% % phaseout.algebraicStates.ePsi_scaled.meas = ePsi;
% % phaseout.algebraicStates.delta_scaled.meas = delta;
