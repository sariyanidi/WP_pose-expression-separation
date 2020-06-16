function [RMSE3D, RMSE2D] = measure_reconstruction_error(S, E, Rs, expression_basis_fname)

T = size(E,1);
if nargin < 3
    Rs = cell(T,1);
    for t=1:T
        Rs{t} = eye(3);
    end
end


X0 = S(1,:)';
Y0 = S(2,:)'; 
Z0 = S(3,:)';

Dc = load(sprintf('./models/%s', expression_basis_fname));

B0 = Dc;


T = size(E,1);
N = length(X0);

RMSE3D = zeros(T, N);
RMSE2D = zeros(T, N);

B = zeros(size(B0));
B(1:3:end,:) = B0(1:N,:);
B(2:3:end,:) = B0(N+1:2*N,:);
B(3:3:end,:) = B0(2*N+1:end,:);





for t=1:T
    Xs_ = B*(E(t,:)');
    Xs = Xs_(1:3:end);
    Ys = Xs_(2:3:end);
    Zs = Xs_(3:3:end);

    X = X0+Xs;
    Y = Y0+Ys;
    Z = Z0+Zs;
    
    %%% Calculate interocular distance
    ix1=37:42; % landmarks corresponding to (camera) left eye
    ix2=43:48; % landmarks corresponding to right eye
    eye1x = mean(X0(ix1));
    eye1y = mean(Y0(ix1));
    eye1z = mean(Z0(ix1));
    
    eye2x = mean(X0(ix2));
    eye2y = mean(Y0(ix2));
    eye2z = mean(Z0(ix2));
    
    iod = norm([eye1x-eye2x; eye1y-eye2y; eye1z-eye2z]);
    
    e = [Xs Ys Zs];
    rmse = sqrt(sum(e.^2,2))/iod;
    RMSE3D(t,:) = rmse;
    
    e = [Xs Ys];
    rmse = sqrt(sum(e.^2,2))/iod;
    RMSE2D(t,:) = rmse;
end





end