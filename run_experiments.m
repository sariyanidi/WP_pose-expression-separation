function run_experiments

close all

addpath('toolbox');
addpath('utils');

% Obtain Fig. 4 and Fig. 7 of paper for pitch rotation
% Figures are stored in folder "results"
[Rrs, Rnrs1, thetas] = one_batch('expbasis_itw.txt', 'pitch');
%[~, Rnrs2] = one_batch('expbasis_basel17.txt', 'pitch');

% Extra (not in CVPR'20 paper): compute errors in pitch angle estimation (in terms of degrees)
clf
for i=1:length(Rrs)
    subplot(1,length(Rrs),i);
    plot(thetas, mean(Rrs{i}), 'b', 'LineWidth', 3);
    hold on
    plot(thetas, mean(Rnrs1{i}), 'r--', 'LineWidth', 3);
    %plot(thetas, mean(Rnrs2{i}), 'k:', 'LineWidth', 3);

    if i==1
        legend({'No PDM', 'Basel-2017 PDM', 'ITWM PDM'});
    end
    pbaspect([1.618 1 1])
    xlim([-45 45]);
end
printpdf(gcf, sprintf('results/rot-errs-pitch-wp.pdf'));


close all

% Obtain Fig. 4 and Fig. 7 of paper for yaw rotation
% Figures are stored in folder "results"
[Rrs1, Rnrs1] = one_batch('expbasis_itw.txt', 'yaw');
%[Rrs2, Rnrs2] = one_batch('expbasis_basel17.txt', 'yaw');

% Extra (not included in CVPR'20 paper): compute errors in yaw angle estimation (in terms of degrees)
clf
for i=1:length(Rrs)
    subplot(1,length(Rrs),i);
    plot(thetas, mean(Rrs{i}), 'b', 'LineWidth', 3);
    hold on
    plot(thetas, mean(Rnrs1{i}), 'r--', 'LineWidth', 3);
    %plot(thetas, mean(Rnrs2{i}), 'k:', 'LineWidth', 3);
    if i==1
        legend({'No PDM', 'Basel-2017 PDM', 'ITWM PDM'});
    end
    
    pbaspect([1.618 1 1])
    xlim([-45 45])
end
printpdf(gcf, sprintf('results/rot-errs-yaw.pdf'));

end





function [ROT_errs_rigid, ROT_errs_nonrigid, thetas] = one_batch(basis_expression_fname, rotation_type)

fprintf('Producing results with expression basis %s\n', basis_expression_fname);
if ~isdir('cache'); mkdir('cache'); end

% ctt: camera transformation type
ctt = 'orthographicxy';

if ~isdir('results'), mkdir('results'); end

% Here we create results for the 6 combinations of 3 fields of view (fovs) and 2 camera distances (tzs)
[~, ~, thetas, mu_RMSE_rigid1, mu_RMSE_expr1, ~, Rr1, Rnr1] = one_fov(1, 1, basis_expression_fname, rotation_type, ctt, 3);
[~, ~, ~, mu_RMSE_rigid2, mu_RMSE_expr2, ~, Rr2, Rnr2] = one_fov(1, 1, basis_expression_fname, rotation_type, ctt, 1.25);
[~, ~, ~, mu_RMSE_rigid3, mu_RMSE_expr3, ~, Rr3, Rnr3] = one_fov(2, 2, basis_expression_fname, rotation_type, ctt, 3);
[~, ~, ~, mu_RMSE_rigid4, mu_RMSE_expr4, ~, Rr4, Rnr4] = one_fov(2, 2, basis_expression_fname, rotation_type, ctt, 1.25);
[~, ~, ~, mu_RMSE_rigid5, mu_RMSE_expr5, ~, Rr5, Rnr5] = one_fov(3, 3, basis_expression_fname, rotation_type, ctt, 3);
[~, ~, ~, mu_RMSE_rigid6, mu_RMSE_expr6, ~, Rr6, Rnr6] = one_fov(3, 3, basis_expression_fname, rotation_type, ctt, 1.25);


% Create and save the panels in Fig. 4 (i.e., 3D-to-2D mapping errors)
close all
subplot(1,6,1);
plot_one_group(mu_RMSE_rigid1, thetas);
legend({'overall', 'brows', 'eyes', 'mouth'});
subplot(1,6,2);
plot_one_group(mu_RMSE_rigid2, thetas);
subplot(1,6,3);
plot_one_group(mu_RMSE_rigid3, thetas);
subplot(1,6,4);
plot_one_group(mu_RMSE_rigid4, thetas);
subplot(1,6,5);
plot_one_group(mu_RMSE_rigid5, thetas);
subplot(1,6,6);
plot_one_group(mu_RMSE_rigid6, thetas);
printpdf(gcf, sprintf('results/Fig4-%s.pdf', rotation_type));


% Create and save the panels in to Fig. 7 (i.e., magnitude spurious expressions)
close all
[~, bn] = fileparts(basis_expression_fname);

subplot(1,6,1);
plot_one_group(mu_RMSE_expr1, thetas);
subplot(1,6,2);
plot_one_group(mu_RMSE_expr2, thetas);
subplot(1,6,3);
plot_one_group(mu_RMSE_expr3, thetas);
subplot(1,6,4);
plot_one_group(mu_RMSE_expr4, thetas);
subplot(1,6,5);
plot_one_group(mu_RMSE_expr5, thetas);
subplot(1,6,6);
plot_one_group(mu_RMSE_expr6, thetas);
printpdf(gcf, sprintf('results/Fig7-%s-%s-wp.pdf', bn, rotation_type));

ROT_errs_rigid = {Rr1,Rr2,Rr3,Rr4,Rr5,Rr6};
ROT_errs_nonrigid = {Rnr1,Rnr2,Rnr3,Rnr4,Rnr5,Rnr6};

end



function plot_one_group(mu_RMSE, thetas)

plot(thetas, mu_RMSE.overall, '-', 'LineWidth', 3);
hold on
plot(thetas, mu_RMSE.brows, ':', 'LineWidth', 3);
plot(thetas, mu_RMSE.eyes, ':.', 'LineWidth', 3);
plot(thetas, mu_RMSE.mouth, '--', 'LineWidth', 3);
xlim([min(thetas), max(thetas)]);

set(gcf,'position',[0,0,750*1.618,750])

ymax = max([max(mu_RMSE.overall),...
    max(mu_RMSE.brows),...
    max(mu_RMSE.eyes),...
    max(mu_RMSE.mouth)]);

plot(zeros(11,1), 0:0.1:1, 'k');

ytick_f = round(ymax, 2);
if ytick_f < ymax
    ytick_f = ymax;
end
ylim([0 ytick_f]);
yticks([0  ytick_f]);

pbaspect([1.618 1 1])

end



function [mmax_rigid, mmax_expr, thetas, mu_RMSE_rigid, mu_RMSE_expr, ...
    fkey, ROT_err_rigid, ROT_err_nonrigid] = one_fov(fovid, tzid, ...
    basis_expression_fname, rotation_type, ctt, dist_coef)

fovs = [30, 60, 90];
tzs = dist_coef*[18, 8.35, 5];

fov = fovs(fovid);
tz  = tzs(tzid);

dist_string = '';
if dist_coef == 1.25
    dist_string = 'close';
else
    dist_string = 'far';
end


fprintf('Producing results for fov=%d degrees (%s distance)\n', fov, dist_string);

fkey = strrep(sprintf('%d-%.2f', fovid, dist_coef), '.', '_');

N = 100;
T = 21;

if nargin < 5
    ctt = 'orthographicxy';
end

cache_fpath = sprintf('cache/exp1_fov=%d-tz=%.1f-basis=%s-rot=%s-ctt=%s-N=%d-T=%d.mat', ...
    fov, tz, basis_expression_fname, rotation_type, ctt, N, T);

if true %  ~exist(cache_fpath)
    % rtt: rigid transform type (it's either rotation or translation)
    rtt = 'rotation';
    
    RMSE_rigid_overall = zeros(N, T);
    RMSE_rigid_eyes = zeros(N, T);
    RMSE_rigid_brows = zeros(N, T);
    RMSE_rigid_mouth = zeros(N, T);
    RMSE_rigid_jaw = zeros(N, T);
    
    RMSE_nonrigid_overall = zeros(N, T);
    RMSE_nonrigid_eyes = zeros(N, T);
    RMSE_nonrigid_brows = zeros(N, T);
    RMSE_nonrigid_mouth = zeros(N, T);
    RMSE_nonrigid_jaw = zeros(N, T);
    
    RMSE_expr_overall = zeros(N, T);
    RMSE_expr_eyes = zeros(N, T);
    RMSE_expr_brows = zeros(N, T);
    RMSE_expr_mouth = zeros(N, T);
    RMSE_expr_jaw = zeros(N, T);
    
    ROT_err_rigid = zeros(N, T);
    ROT_err_nonrigid = zeros(N, T);
    
    [~,~,~,~,~,~,~,~,thetas] = measure_RMSE_error_on_rotated_sequence(1, ...
        fov, tz, T, rtt, ctt, rotation_type,  basis_expression_fname);
    
    tic
    for n=1:N
        if mod(n,20) == 0
            fprintf('\t Processed %03d/%d faces (in %.2f secs) \n', n, N, toc);
        end
        [~, ~, rmse_rigid, rmse_nonrigid, ~,rmse_expr3D, rot_err_rigid, rot_err_nonrigid, ~] = ...
            measure_RMSE_error_on_rotated_sequence(n, fov, tz, T, rtt, ctt, rotation_type,  basis_expression_fname);
        RMSE_rigid_overall(n,:) = rmse_rigid.overall;
        RMSE_rigid_eyes(n,:) = rmse_rigid.eyes;
        RMSE_rigid_brows(n,:) = rmse_rigid.brows;
        RMSE_rigid_mouth(n,:) = rmse_rigid.mouth;
        RMSE_rigid_jaw(n,:) = rmse_rigid.jaw;
        
        RMSE_nonrigid_overall(n,:) = rmse_nonrigid.overall;
        RMSE_nonrigid_eyes(n,:) = rmse_nonrigid.eyes;
        RMSE_nonrigid_brows(n,:) = rmse_nonrigid.brows;
        RMSE_nonrigid_mouth(n,:) = rmse_nonrigid.mouth;
        RMSE_nonrigid_jaw(n,:) = rmse_nonrigid.jaw;
        
        RMSE_expr_overall(n,:) = rmse_expr3D.overall;
        RMSE_expr_eyes(n,:) = rmse_expr3D.eyes;
        RMSE_expr_brows(n,:) = rmse_expr3D.brows;
        RMSE_expr_mouth(n,:) = rmse_expr3D.mouth;
        RMSE_expr_jaw(n,:) = rmse_expr3D.jaw;

        
        ROT_err_rigid(n,:) = rot_err_rigid;
        ROT_err_nonrigid(n,:) = rot_err_nonrigid;
    end    
    save(cache_fpath);
else
    load(cache_fpath);
end

mu_RMSE_rigid.overall = mean(RMSE_rigid_overall);
mu_RMSE_rigid.brows = mean(RMSE_rigid_brows);
mu_RMSE_rigid.eyes = mean(RMSE_rigid_eyes);
mu_RMSE_rigid.mouth = mean(RMSE_rigid_mouth);

mu_RMSE_expr.overall = mean(RMSE_expr_overall);
mu_RMSE_expr.brows = mean(RMSE_expr_brows);
mu_RMSE_expr.eyes = mean(RMSE_expr_eyes);
mu_RMSE_expr.mouth = mean(RMSE_expr_mouth);

mmax_rigid = max([max(mu_RMSE_rigid.overall),...
    max(mu_RMSE_rigid.brows),...
    max(mu_RMSE_rigid.eyes),...
    max(mu_RMSE_rigid.mouth)]);

mmax_expr = max([max(mu_RMSE_expr.overall),...
    max(mu_RMSE_expr.brows),...
    max(mu_RMSE_expr.eyes),...
    max(mu_RMSE_expr.mouth)]);

pbaspect([1.618 1 1])

end
