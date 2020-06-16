function [q, tr, e, rmse] = fit_3d_model_withexpressions_gaussnewton(camera, S, x, y, q00, q10, q20, q30, tr0, e0, expresion_basis_fname, frame_id)

if nargin < 12
    frame_id = 0;
end
% First run framewise_head_movement; this will:
% Take X, Y, Z from Tomasi-Kanade
% Initialize q0, q1, q2, q3 by solving PnP 

% load('D1.mat');
% B = Dc;

% phi = 9999;

warning('off', 'MATLAB:nearlySingularMatrix');

xnorm = x-mean(x);
ynorm = -(y-mean(y));
N = length(xnorm);

% % % % % % % Dc = load('expression_basis_c.txt');
Dc = load(sprintf('./models/%s', expresion_basis_fname));
B0 = Dc;
B = zeros(size(B0));
B(1:3:end,:) = B0(1:N,:);
B(2:3:end,:) = B0(N+1:2*N,:);
B(3:3:end,:) = B0(2*N+1:end,:);
% B = B0;
e0 = zeros(size(B,2),1);


X = S(1,:);
Y = S(2,:); 
Z = S(3,:); 

cx = 250;
cy = 250;


%%

q0 = q00;
q1 = q10;
q2 = q20;
q3 = q30;
tr = tr0;
e = e0;



% q0 = 1;
% q1 = 0;
% q2 = 0;
% q3 = 0;

sigmax = camera.sigmax;
sigmay = camera.sigmay;

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
    
    RN = cell(N,1);
    dRq1N = cell(N,1);
    dRq2N = cell(N,1);
    dRq3N = cell(N,1);
    for n=1:N
        RN{n} = R(1:2,:);
        RNfull{n} = R;
        dRq1N{n} = dRq1;
        dRq2N{n} = dRq2;
        dRq3N{n} = dRq3;
    end
    
    RB = blkdiag(RN{:});
    dR1 = blkdiag(dRq1N{:});
    dR2 = blkdiag(dRq2N{:});
    dR3 = blkdiag(dRq3N{:});
    Rfull = blkdiag(RNfull{:});
    
    P = [X; Y; Z];
    P = P(:);
    
    sigmaP = [sigmax*X; sigmay*Y; sigmax*Z];
    sigmaP = sigmaP(:);
    
    
    
    Ps = B*e;
    
    sigmaPs = Ps;
    sigmaPs(1:3:end) = sigmax*sigmaPs(1:3:end); 
    sigmaPs(2:3:end) = sigmay*sigmaPs(2:3:end); 
    sigmaPs(3:3:end) = sigmay*sigmaPs(3:3:end); 
    
    dP_q1 = dR1*sigmaP+dR1*sigmaPs;
    dP_q2 = dR2*sigmaP+dR2*sigmaPs;
    dP_q3 = dR3*sigmaP+dR3*sigmaPs;
    
    v = Rfull*P + Rfull*Ps + repmat(tr, [N,1]);
    vX = v(1:3:end)';
    vY = v(2:3:end)';
    vZ = v(3:3:end)';

    dP_dtX =  sigmax*[ones(1, N); zeros(1, N)];
    dP_dtY =  sigmay*[zeros(1, N); ones(1, N)];
    
    A = RB*B;
    A(1:2:end) = sigmax*A(1:2:end);
    A(2:2:end) = sigmay*A(2:2:end); 

    
    dPx_dsx = vX;
    dPy_dsx = 0*vX;
    dPx_dsy = 0*vY;
    dPy_dsy = vY;
    dP_dsx = [dPx_dsx; dPy_dsx];
    dP_dsy = [dPx_dsy; dPy_dsy];
    
    J = [dP_dsx(:), dP_dsy(:), dP_q1, dP_q2, dP_q3, dP_dtX(:), dP_dtY(:), A];
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
    
    if (norm(delta_c) < 0.001) || (numiter == 10)
        break
    end
    
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    
    tr(1:2) = tr(1:2) + delta_c(6:7);
    delta_e = delta_c(8:8+length(e)-1);
    e = e+delta_e;
    
    sigmax = sigmax + delta_c(1);
    sigmay = sigmay + delta_c(2);
end


err = [xnorm; ynorm] - [sigmax.*vX; sigmay.*vY;];
rmse = sqrt(sum(err.^2));

if camera.store == true
    if ~exist('./figures')
        mkdir('./figures');
    end
    [~, bn] = fileparts(expresion_basis_fname);
    fpath = sprintf('figures/%02d-nonrigid-%s-%s.pdf', frame_id, bn);
    
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
nose2 = [32, 33, 34, 35, 36];
 
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




