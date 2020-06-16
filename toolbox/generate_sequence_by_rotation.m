function [pX, pY, Xo, Yo, Zo, thetas, iod] = generate_sequence_by_rotation(S, T, camera, rotation_type)


N = size(S,2);

X = [S(1,:); S(2,:); S(3,:)];
pX = zeros(N, T);
pY = zeros(N, T);
Xo = zeros(N, T);
Yo = zeros(N, T);
Zo = zeros(N, T);

P = X(:);

tr = repmat([0 0 -camera.tz]', N, 1);
theta_min = -45;
theta_max = 45;
thetas = theta_min:((theta_max-theta_min)/(T-1)):theta_max;

iod = 0; % interocular distance
for t=1:T
    if strcmp(rotation_type, 'yaw')
        R = eul2rotm([0 deg2rad(thetas(t)) 0]);
    else
        R = eul2rotm([0 0 deg2rad(thetas(t))]);
    end
    RN = cell(N,1);
    Is = cell(N,1);
    for n=1:N
        RN{n} = R;
        Is{n} = (camera.Intrinsic)';
    end
    RN = blkdiag(RN{:});
    CI = blkdiag(Is{:});
    V = RN*P+tr;
    xp = CI*V;
    x = xp(1:3:end);
    y = xp(2:3:end);
    z = xp(3:3:end);
    
    x = x./z;
    y = y./z;

    vx = V(1:3:end);
    vy = V(2:3:end);

    vz = V(3:3:end);
    Xo(:,t) = vx;
    Yo(:,t) = vy;
    Zo(:,t) = vz;
    

    if strcmp(camera.type, 'perspective')
        pX(:,t) = x;
        pY(:,t) = -y;
        
    elseif strcmp(camera.type, 'orthographic')
        pX(:,t) = camera.sigma*vx;
        pY(:,t) = -(camera.sigma*vy);
    end
    
    if camera.visualize
        plot(x,y, 'x');
        axis square

    end
    
    %%% Calculate interocular distance
    ix1=37:42; % landmarks corresponding to (camera) left eye
    ix2=43:48; % landmarks corresponding to right eye
    eye1x = mean(x(ix1)); eye1y = mean(y(ix1));
    eye2x = mean(x(ix2)); eye2y = mean(y(ix2));
    iod = norm([eye1x-eye2x; eye1y-eye2y]);
end

pX = pX';
pY = pY';



end