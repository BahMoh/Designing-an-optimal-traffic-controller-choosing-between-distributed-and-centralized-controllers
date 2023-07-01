clc;
clear;
close all;

A1 = [-3 -2;1 0];
A2 = [-5.5 -3;2 0];
A3 = [-4 -1.75;1 0];
A4 = [-4 -2;2 0];

B1 = [4;0];
B2 = [2;0];
B3 = [1;0];
B4 = [2;0];

R1 = 1;
R2 = 0.5;
R3 = 1;
R4 = 2;

Q1 = [1 0;0 1];
Q2 = [1 0;0 2];
Q3 = [2 0;0 1];
Q4 = [2 0;0 2];

Q = [Q1 Q2 Q3 Q4];
R = [R1 R2 R3 R4];

k1s = dlqr(A1,B1,Q1,R1);
k2s = dlqr(A2,B2,Q2,R2);
k3s = dlqr(A3,B3,Q3,R3);
k4s = dlqr(A4,B4,Q4,R4);

k1w = dlqr(A1,B1,Q4,R4);
k2w = dlqr(A2,B2,Q1,R1);
k3w = dlqr(A3,B3,Q2,R2);
k4w = dlqr(A4,B4,Q3,R3);

samples = 10;

[u1s,x1s] = controlsignalState(A1,B1,k1s,samples);
[u2s,x2s] = controlsignalState(A2,B2,k2s,samples);
[u3s,x3s] = controlsignalState(A3,B3,k3s,samples);
[u4s,x4s] = controlsignalState(A4,B4,k4s,samples);

[u1w,x1w] = controlsignalState(A1,B1,k1w,samples);
[u2w,x2w] = controlsignalState(A2,B2,k2w,samples);
[u3w,x3w] = controlsignalState(A3,B3,k3w,samples);
[u4w,x4w] = controlsignalState(A4,B4,k4w,samples);
%% Allcombs
U1 = [u1s;u1w];
U2 = [u2s;u2w];
U3 = [u3s;u3w];
U4 = [u4s;u4w];

X1 = [x1s;x1w];
X2 = [x2s;x2w];
X3 = [x3s;x3w];
X4 = [x4s;x4w];

combinations = allcomb([1     2],[1     2],[1     2],[1     2]);

U_1 = [];
U_2 = [];
U_3 = [];
U_4 = [];

J = 0;
JJ = zeros(length(combinations),1);

for k = 1:length(combinations)
    U_1 = U1(combinations(k,1),:);
    U_2 = U2(combinations(k,2),:);
    U_3 = U3(combinations(k,3),:);
    U_4 = U4(combinations(k,4),:);
    
    X_1 = X1(2*combinations(k,1)-1:2*combinations(k,1),:);
    X_2 = X2(2*combinations(k,2)-1:2*combinations(k,2),:);
    X_3 = X3(2*combinations(k,3)-1:2*combinations(k,3),:);
    X_4 = X4(2*combinations(k,4)-1:2*combinations(k,4),:);

    for j = 1:length(u1s)
        Sum = (X_1(:,j)'*Q1*X_1(:,j) + U_1(:,j)'*R1*U_1(:,j))+...
              (X_2(:,j)'*Q2*X_2(:,j) + U_2(:,j)'*R2*U_2(:,j))+...
              (X_3(:,j)'*Q3*X_3(:,j) + U_3(:,j)'*R3*U_3(:,j))+...
              (X_4(:,j)'*Q4*X_4(:,j) + U_4(:,j)'*R4*U_4(:,j));
        J = Sum + J;
    end
    JJ(k,1) = J;
end
[minimum, index] = min(JJ);

u1 = U1(combinations(index,1), :);
u2 = U2(combinations(index,1), :);
u3 = U3(combinations(index,1), :);
u4 = U4(combinations(index,1), :);

x1 = X1(2*combinations(index,1)-1:2*combinations(index,1),:);
x2 = X2(2*combinations(index,1)-1:2*combinations(index,1),:);
x3 = X3(2*combinations(index,1)-1:2*combinations(index,1),:);
x4 = X4(2*combinations(index,1)-1:2*combinations(index,1),:);

%% Plot the results
figure(1)
subplot(2,2,1);
plot(x1(1,:))
hold on
plot(x1(2,:))
legend('sub1 state 1', 'sub1 state 2');
title("subsystem 1");
xlabel("samples")
ylabel("states for subsystem1")

subplot(2,2,2);
plot(x2(1,:))
hold on
plot(x2(2,:))
legend('sub2 state 1', 'sub2 state 2');
title("subsystem 2");
xlabel("samples")
ylabel("states for subsystem2")

subplot(2,2,3);
plot(x3(1,:))
hold on
plot(x3(2,:))
legend('sub3 state 1', 'sub3 state 2');
title("subsystem 3");
xlabel("samples")
ylabel("states for subsystem3")

subplot(2,2,4);
plot(x4(1,:))
hold on
plot(x4(2,:))
legend('sub4 state 1', 'sub4 state 2');
title("subsystem 4");
xlabel("samples")
ylabel("states for subsystem4")

%Plot control signals
figure(2)
subplot(2,2,1);
plot(u1)
legend("sub 1 control signal")
title("subsystem 1")
xlabel("samples")
ylabel("control signal for subsystem 1")

subplot(2,2,2);
plot(u2)
legend("sub 2 control signal")
title("subsystem 2")
xlabel("samples")
ylabel("control signal for subsystem 2")

subplot(2,2,3);
plot(u3)
legend("sub 3 control signal")
title("subsystem 3")
xlabel("samples")
ylabel("control signal for subsystem 3")

subplot(2,2,4);
plot(u4)
legend("sub 4 control signal")
title("subsystem 4")
xlabel("samples")
ylabel("control signal for subsystem 4")