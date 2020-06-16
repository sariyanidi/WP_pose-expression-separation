function [pose_rigid, pose_nonrigid, E, sigmaX, sigmaY, errs_rigid, errs_nonrigid, Rs] = ...
    compute_coefficients_and_errors(S, pX, pY, camera, transform_type, expresion_basis_fname, iod)

% % % % % % Dc = load('expression_basis_c.txt');
Dc = load(sprintf('./models/%s', expresion_basis_fname));
B = Dc;
% e0 = zeros(size(B,2),1);
K = size(B,2);
T = size(pX,1);
E = zeros(T, K);
sigmaX = zeros(T,1);
sigmaY = zeros(T,1);

yaw_rigid = zeros(T,1);
pitch_rigid = zeros(T,1);
roll_rigid = zeros(T,1);

yaw_nonrigid = zeros(T,1);
pitch_nonrigid = zeros(T,1);
roll_nonrigid = zeros(T,1);

err_rigid_overall = zeros(T,1);
err_rigid_eyes = zeros(T,1);
err_rigid_brows = zeros(T,1);
err_rigid_mouth = zeros(T,1);
err_rigid_jaw = zeros(T,1);

err_nonrigid_overall = zeros(T,1);
err_nonrigid_eyes = zeros(T,1);
err_nonrigid_brows = zeros(T,1);
err_nonrigid_mouth = zeros(T,1);
err_nonrigid_jaw = zeros(T,1);

Rs = cell(T,1);


for t=1:T
    x = pX(t,:);
    y = pY(t,:);
    
    [q0, tr0, camera, rmse_errs] = fit_3d_model_gaussnewton(camera, S, x, y, 1, 0, 0, 0, [0;0;1], t);
    [q, tr, e, rmse_errs_exp] = fit_3d_model_withexpressions_gaussnewton(camera, S, x, y, q0(1), q0(2), q0(3), q0(4), tr0, zeros(3,1), expresion_basis_fname, t);
    
    sigmaX(t) = camera.sigmax;
    sigmaY(t) = camera.sigmay;
    
    % indices for the following facial features
    eyes = 37:48;
    brows = 18:27;
    mouth = 49:68;
    jaw = 1:17;
    
    err_rigid_overall(t) = mean(rmse_errs);
    err_rigid_eyes(t) = mean(rmse_errs(eyes));
    err_rigid_brows(t) = mean(rmse_errs(brows));
    err_rigid_mouth(t) = mean(rmse_errs(mouth));
    err_rigid_jaw(t) = mean(rmse_errs(jaw));
    
    err_nonrigid_overall(t) = mean(rmse_errs_exp);
    err_nonrigid_eyes(t) = mean(rmse_errs_exp(eyes));
    err_nonrigid_brows(t) = mean(rmse_errs_exp(brows));
    err_nonrigid_mouth(t) = mean(rmse_errs_exp(mouth));
    err_nonrigid_jaw(t) = mean(rmse_errs_exp(jaw));
    
    E(t,:) = e;
    
    eul_rigid = rotm2eul(quat2rotm(q0));
    pitch_rigid(t) = eul_rigid(3);
    yaw_rigid(t) = eul_rigid(2);
    
    eul_nonrigid = rotm2eul(quat2rotm(q));
    pitch_nonrigid(t) = eul_nonrigid(3);
    yaw_nonrigid(t) = eul_nonrigid(2);
end




pose_rigid = [yaw_rigid pitch_rigid roll_rigid];
pose_nonrigid = [yaw_nonrigid pitch_nonrigid roll_nonrigid];


errs_rigid.overall = err_rigid_overall/iod;
errs_rigid.eyes = err_rigid_eyes/iod;
errs_rigid.brows = err_rigid_brows/iod;
errs_rigid.mouth = err_rigid_mouth/iod;
errs_rigid.jaw = err_rigid_jaw/iod;


errs_nonrigid.overall = err_nonrigid_overall/iod;
errs_nonrigid.eyes = err_nonrigid_eyes/iod;
errs_nonrigid.brows = err_nonrigid_brows/iod;
errs_nonrigid.mouth = err_nonrigid_mouth/iod;
errs_nonrigid.jaw = err_nonrigid_jaw/iod;

