function [r] = quat2rotm(q)
    % This function is taken from https://github.com/robotology/mex-wholebodymodel
    % In fact the function comes by default in not-so-old versions of MATLAB but we include it here for completeness

    %Check if vector and if have it has enough input values
    if ~isvector(q) || (length(q))~=4 
        error('Input to quat2r must be a vector and must have 4 values')
    end
    
    %Initialize
    r=zeros(3,3);
    
    %Normailize Quaternion before using
    q=q/norm(q);
    
    %Quaternion components
    qw=q(1);
    qx=q(2);
    qy=q(3);
    qz=q(4);
    
    %Rotation Matrix Elements
    r(1,1) = 1 - 2*qy^2 - 2*qz^2;
    r(1,2) = 2*qx*qy - 2*qz*qw;
    r(1,3) = 2*qx*qz + 2*qy*qw;
    
    r(2,1) = 2*qx*qy + 2*qz*qw;
    r(2,2) = 1 - 2*qx^2 - 2*qz^2;
    r(2,3) = 2*qy*qz - 2*qx*qw;
    
    r(3,1) = 2*qx*qz - 2*qy*qw;
    r(3,2) = 2*qy*qz + 2*qx*qw;
    r(3,3) = 1 - 2*qx^2 - 2*qy^2;
    
    %Equivalent built-in MATLAB function: quat2rotm
end
