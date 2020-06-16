function [S, E, rmse_rigid, rmse_nonrigid, rmse_expr2D, rmse_expr3D, rot_err_rigid, rot_err_nonrigid, thetas] = ...
    measure_RMSE_error_on_rotated_sequence(subject_id, fov, tz, T, rigid_transform_type, camera_transform_type, rotation_type, expression_basis_fname)

%%% Step 1) START: Generate the rotated sequences of points
load(sprintf('facial_shapes/S_%03d.mat', subject_id));
S = S';

imw = 640;
imh = 480;

S(1,:) = S(1,:)-mean(S(1,:));
S(2,:) = S(2,:)-mean(S(2,:));
S(3,:) = S(3,:)-mean(S(3,:));
S = S/std(S(:));

camera.type = 'perspective';
camera.phi = 999;
camera.tz = tz;

camera.visualize = false;
camera.store = false;
% if subject_id == 1
%     camera.store = true;  
% end
% parameters of the camera, needed to generate the true points
phix = (imw/2)/tand(fov/2);
phiy = (imh/2)/tand(fov/2);

camera.Intrinsic = eye(3);
camera.Intrinsic(3,2) = imh/2;
camera.Intrinsic(3,1) = imw/2;
camera.Intrinsic(1,1) = phix;
camera.Intrinsic(2,2) = phiy;
camera.sigma = phix/tz;


outdir = 'coeffs';

if strcmp(rigid_transform_type, 'rotation')
    seq_name = sprintf('seq-S%03d-%d-%.1f-%s-%s-T%02d', subject_id, fov, tz, camera_transform_type, rotation_type, T);
elseif strcmp(rigid_transform_type, 'translation')
    seq_name = sprintf('seq-S%03d-%d-%.1f-%s-%s-T%02d', subject_id, fov, tz, camera_transform_type, direction, T);
end

seq_path = sprintf('%s/%s', outdir, seq_name);

if ~isdir(outdir), mkdir(outdir); end
if ~isdir(seq_path), mkdir(seq_path); end

% Generate the rotated sequence of points
if strcmp(rigid_transform_type, 'rotation')
    [pX, pY, vX, vY, vZ, thetas, iod] = generate_sequence_by_rotation(S, T, camera, rotation_type);
elseif strcmp(rigid_transform_type, 'translation')
    % to be implemented
end
%%% Step 1) END: generate the rotated sequences of points



%%% Step 2) START: fit 3d morphable models (PDMs) and also compute 3D-to-2D mapping errors
transform_type = camera_transform_type;

[pose_rigid, pose_nonrigid, E, ~, ~, rmse_rigid, rmse_nonrigid, Rs] = ...
    compute_coefficients_and_errors(S, pX, pY, camera, transform_type, expression_basis_fname, iod);

 
angle_id = 1;
if strcmp(rotation_type, 'pitch')
    angle_id = 2;
end

rot_err_rigid = abs(rad2deg(pose_rigid(:,angle_id))-thetas');
rot_err_nonrigid = abs(rad2deg(pose_nonrigid(:,angle_id))-thetas');

for t=1:T
    e = E(t,:)';
    
    fpath = sprintf('%s/%04d_%s-shape.txt', seq_path,  t, strrep(expression_basis_fname, '.txt', ''));
    save(fpath, 'e', '-ASCII')
end
%%% Step 2) END: fit 3d morphable models (PDMs) and also compute 3D-to-2D mapping errors

%%% Step 3) START: Compute magnitude of spurious expressions
[RMSE3D, RMSE2D] = measure_reconstruction_error(S, E, Rs, expression_basis_fname);

eyesix = 37:48;
browsix = 18:27;
mouthix = 49:68;
jawix = 1:17;

rmse_expr3D.overall = mean(RMSE3D');
rmse_expr3D.brows = mean(RMSE3D(:,browsix)');
rmse_expr3D.eyes = mean(RMSE3D(:,eyesix)');
rmse_expr3D.mouth = mean(RMSE3D(:,mouthix)');
rmse_expr3D.jaw = mean(RMSE3D(:,jawix)');

rmse_expr2D.overall = mean(RMSE2D');
rmse_expr2D.brows = mean(RMSE2D(:,browsix)');
rmse_expr2D.eyes = mean(RMSE2D(:,eyesix)');
rmse_expr2D.mouth = mean(RMSE2D(:,mouthix)');
rmse_expr2D.jaw = mean(RMSE2D(:,jawix)');

%%% Step 3) END: Compute magnitude of spurious expressions

