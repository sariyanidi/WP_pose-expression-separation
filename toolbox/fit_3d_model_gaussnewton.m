function [q, tr, camera, rmse] = fit_3d_model_gaussnewton(camera, S, x, y, q00, q10, q20, q30, tr0, frame_id)

if nargin < 10
    frame_id = 0;
end

% First run framewise_head_movement; this will:
% Take X, Y, Z from Tomasi-Kanade

% Initialize q0, q1, q2, q3 by solving PnP 

% phi = 9999;
S = double(S);
X = S(1,:);
Y = S(2,:); 
Z = S(3,:);

xnorm = x-mean(x);
ynorm = -(y-mean(y));
N = length(xnorm);

%%

q0 = q00;
q1 = q10;
q2 = q20;
q3 = q30;
tr = tr0;

sigmax = camera.sigma;
sigmay = camera.sigma;
sigma = camera.sigma;

numiter = 0;

while true
    
    numiter = numiter+1;

    R = quat2rotm([q0, q1, q2, q3]);
    dRq1 = [    0,  2*q2,  2*q3; 2*q2, -4*q1, -2*q0; 2*q3,  2*q0, -4*q1];
    dRq2 = [ -4*q2, 2*q1,  2*q0; 2*q1,    0,  2*q3; -2*q0, 2*q3, -4*q2]; 
    dRq3 = [ -4*q3, -2*q0, 2*q1; 2*q0, -4*q3, 2*q2; 2*q1,  2*q2,    0];
    dRq1 = dRq1(1:2,:);
    dRq2 = dRq2(1:2,:);
    dRq3 = dRq3(1:2,:);
    
    dRq1N = cell(N,1);
    dRq2N = cell(N,1);
    dRq3N = cell(N,1);
    for n=1:N
        dRq1N{n} = dRq1;
        dRq2N{n} = dRq2;
        dRq3N{n} = dRq3;
    end
    
    dR1 = blkdiag(dRq1N{:});
    dR2 = blkdiag(dRq2N{:});
    dR3 = blkdiag(dRq3N{:});
    
    P = [X; Y; Z];
    sigmaP = [sigmax*X; sigmay*Y; sigmax*Z];
    sigmaP = sigmaP(:);
    P = P(:);
    
    dP_q1 = dR1*sigmaP;
    dP_q2 = dR2*sigmaP;
    dP_q3 = dR3*sigmaP;
    
    v = R*[X; Y; Z]+repmat(tr,[1, N]);
    vX = v(1,:);
    vY = v(2,:);
    vZ = v(3,:);

    dP_dtX =  sigmax*[ones(1, N); zeros(1, N)];
    dP_dtY =  sigmay*[zeros(1, N); ones(1, N)];
    
    dPx_dsx = vX;
    dPy_dsx = 0*vX;
    dPx_dsy = 0*vY;
    dPy_dsy = vY;
    dP_dsx = [dPx_dsx; dPy_dsx];
    dP_dsy = [dPx_dsy; dPy_dsy];
    
    J = [dP_dsx(:), dP_dsy(:), dP_q1, dP_q2, dP_q3, dP_dtX(:), dP_dtY(:)];
    A = J'*J;
    
    
    if camera.visualize
        clf
        plot_single(xnorm, ynorm, 'b');
        hold on
        plot_single(sigmax.*vX, sigmay.*vY, 'r--');
    end

    b = [xnorm; ynorm] - [sigmax.*vX; sigmay.*vY;];
    b = b(:);
    b = J'*b;
    delta_c = A\b;
    delta_q13 = delta_c(3:5);
    delta_q = [1; delta_q13];
    delta_q = delta_q/norm(delta_q);


    q = quatmultiply(delta_q', [q0, q1, q2, q3]);
    qnorm = norm(q-[q0, q1, q2, q3]);
    
    if (norm(delta_c) < 0.01) || (numiter == 50)
        break
    end
    
    
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    sigmax = sigmax + delta_c(1);
    sigmay = sigmay + delta_c(2);
    tr(1:2) = tr(1:2) + delta_c(6:end);
end

camera.sigma = sigma;
camera.sigmax = sigmax;
camera.sigmay = sigmay;

err = [xnorm; ynorm] - [sigmax.*vX; sigmay.*vY;];
rmse = sqrt(sum(err.^2));

if camera.store == true
    if ~exist('./figures')
        mkdir('./figures');
    end
    fpath = sprintf('figures/%02d.pdf', frame_id);
    
    if ~exist(fpath)
        printpdf(gcf, fpath);
    end
end



function plot_single(x, y, plotopts)
y = -y;
% clear
 
 
% close all
 
pcolor = [0, 0.4470, 0.7410];
bcolor = [0.8500, 0.3250, 0.0980]; % brows color
ecolor = [0.9290, 0.6940, 0.1250]; % eyes color
mcolor = [0.4940, 0.1840, 0.5560]; % mouth color
rcolor = [0.4660, 0.6740, 0.1880]; % color of the rest
 
% plot(S(:,1),S(:,2), 'Color', pcolor);
hold on
leye = [37:42 37];
reye = [43:48 43];
lbrow = 18:22;
rbrow = 23:27;
mouth = 49:68;
jaw = 1:17;
 
lw = 2;
 
ms = 30;
% plot(S(:,1),S(:,2), 'x');
hold on
plot(x(leye), y(leye), plotopts, 'LineWidth', lw);
% plot(x(leye), y(leye), '.', 'Color', ecolor, 'LineWidth', lw, 'MarkerSize', ms);
 
plot(x(reye), y(reye), plotopts, 'LineWidth', lw);
% plot(x(reye), y(reye),  '.', 'Color', ecolor, 'LineWidth', lw, 'MarkerSize', ms);
 
plot(x(lbrow), y(lbrow), plotopts, 'LineWidth', lw);
% plot(x(lbrow), y(lbrow),  '.', 'Color', bcolor, 'LineWidth', lw, 'MarkerSize', ms);
 
plot(x(rbrow), y(rbrow), plotopts, 'LineWidth', lw);
% plot(x(rbrow), y(rbrow),  '.', 'Color', bcolor, 'LineWidth', lw, 'MarkerSize', ms);
 
plot(x(jaw), y(jaw), plotopts, 'LineWidth', lw);
% plot(x(jaw), y(jaw),  '.', 'Color', rcolor, 'LineWidth', lw, 'MarkerSize', ms);
 
plot(x(mouth), y(mouth), plotopts, 'LineWidth', lw);
% plot(x(mouth), y(mouth),  '.', 'Color', mcolor, 'LineWidth', lw, 'MarkerSize', ms);
 
nose1 = [28:31 34];
nose2 = [32,33,34,35,36];
 
plot(x(nose1), y(nose1), plotopts, 'LineWidth', lw);
% plot(x(nose1), y(nose1),  '.', 'Color', rcolor, 'LineWidth', lw, 'MarkerSize', ms);
 
plot(x(nose2), y(nose2), plotopts, 'LineWidth', lw);
% plot(x(nose2), y(nose2),  '.', 'Color', rcolor, 'LineWidth', lw, 'MarkerSize', ms);

xlim([-100 100]);
ylim([-100 100]);
axis equal
pbaspect([1 1 1])
axis off

% axis equal
% axis square
% % axis off




