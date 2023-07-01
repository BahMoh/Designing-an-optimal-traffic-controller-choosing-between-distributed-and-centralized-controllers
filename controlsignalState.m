% function [u,x] = controlsignalState(A,B,k,samples)
% %CONTROLSIGNALSTATE Summary of this function goes here
% %   Detailed explanation goes here
% x0 = [1.5;1.5];
% x = zeros(2,10);
% x(:,1) = x0;
% u = zeros(1,10);
% u(1,1) = 0;
% for i = 1:samples - 1
%     x(:,i+1) = A*x(:,i) + B*(-k)*x(:,i);
%     u(i+1) = -k*x(:,i);
% end
% end
function [u,x] = controlsignalState(A,B,k,samples)
%CONTROLSIGNALSTATE Summary of this function goes here
%   Detailed explanation goes here
x0 = [1.5;1.5];
x = zeros(2,10);
x(:,1) = x0;
u = zeros(1,10);
u(1,1) = -k*x(:,1);
for i = 1:samples - 1
    x(:,i+1) = A*x(:,i) + B*(-k)*x(:,i);
    u(i+1) = -k*x(:,i+1);
end
end

