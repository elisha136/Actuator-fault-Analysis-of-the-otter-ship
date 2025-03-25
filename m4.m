%% Otter Actuator Fault Diagnosis using AEKF (Modified Solution)
% This script implements an adaptive extended Kalman filter (AEKF) for diagnosing
% actuator faults on a vessel with enhanced control scheduling and improved visualization.

clear; clc;

%% Parameters
Xu     = -0.7225;   Xuu    = -1.3274;   m      = 23.8;     X_udot = -2.0;
Yv     = -0.8612;   Yvv    = -36.2823;  Yv_dot = -10.0;    xg     = 0.046;
Yr     = 0.1079;    Nrr    = -1;        Yr_dot = 0.0;      Nv_dot = 0.0;
Nv     = 0.1052;    Nr     = -0.5;      Izz    = 1.76;     Nr_dot = -1.0;

%% Time Variables
dt = 0.001;
T  = 15;
t  = 0:dt:T;    % Start at 0 to allow valid indexing

%% Noise Matrices
QF = 10 * eye(6);
RF = 40 * eye(6);

%% Initial Values
x   = 0; y   = 0; yaw = 0;
u   = 2; v   = 1; r   = 0;
tau_u = -30; tau_v = 0; tau_r = 0;
eta   = [x; y; yaw];            % vessel position and heading
nu    = [u; v; r];              % body-fixed velocities
tau_c = [tau_u; tau_v; tau_r];   % initial control input
tau_true = tau_c;
X     = [eta; nu];              % full state vector
X_true = X;
xhat    = zeros(6,1);           % initial state estimate
theta   = zeros(3,1);           % true fault parameter
thetahat = zeros(3,1);          % fault parameter estimate

%% Modified Control Unit with Smooth Transitions
% Define piecewise control inputs with smooth transitions using cosine blending
t_control = [0, 3, 6, 9, 12, 15];             % Control change timepoints (s)
tau_u_profile = [ -30, 200, 50, 100, -50, 0 ];  % Surge control profile
tau_v_profile = [ 0, 50, -40, 30, -20, 0 ];     % Sway control profile
tau_r_profile = [ 0, 0.2, -0.5, 0.3, -0.2, 0 ];  % Yaw control profile

% Preallocate control vectors
tau_u = zeros(size(t));
tau_v = zeros(size(t));
tau_r = zeros(size(t));

% Generate smooth control signals using cosine interpolation
for i = 1:length(t_control)-1
    idx = t >= t_control(i) & t < t_control(i+1);
    tau_u(idx) = (tau_u_profile(i+1) - tau_u_profile(i))/2 * ...
        (1 - cos(pi*(t(idx)-t_control(i))/(t_control(i+1)-t_control(i)))) + tau_u_profile(i);
    tau_v(idx) = (tau_v_profile(i+1) - tau_v_profile(i))/2 * ...
        (1 - cos(pi*(t(idx)-t_control(i))/(t_control(i+1)-t_control(i)))) + tau_v_profile(i);
    tau_r(idx) = (tau_r_profile(i+1) - tau_r_profile(i))/2 * ...
        (1 - cos(pi*(t(idx)-t_control(i))/(t_control(i+1)-t_control(i)))) + tau_r_profile(i);
end

tau_c = [tau_u; tau_v; tau_r]; % Full control input vector based on smooth transitions

%% Tuning Parameters
lambda = 0.995;
a      = 0.999;

%% System Matrices
A = eye(6);
M = [ m - X_udot,         0,           0;
      0,           m - Yv_dot,  m*xg - Yr_dot;
      0,           m*xg - Nv_dot, Izz - Nr_dot ];
m11 = M(1,1); m22 = M(2,2); m23 = M(2,3); m32 = M(3,2);
B_c = [ zeros(3); inv(M)*eye(3) ];  % control input matrix
C   = eye(6);
Psi_d = -dt * B_c * diag(tau_c);     % fault injection mapping

%% Estimation Parameters
Pplus = eye(rank(A));
S = 0.1 * eye(3);
UpsilonPlus = 0 * B_c;  % zero matrix (same size as B_c)

%% Preallocate Arrays for Plotting
x_true_vec   = zeros(6, length(t));
eta_vec      = zeros(3, length(t));
nu_vec       = zeros(3, length(t));
tau_c_vec    = zeros(3, length(t));
theta_vec    = zeros(3, length(t));
thetahat_vec = zeros(3, length(t));
xhat_vec     = zeros(6, length(t));

%% True State Simulation
for i = 1:length(t)
    % Schedule a change in true actuation at 25% and 50% of the simulation:
    if i == round(0.25 * length(t))
        tau_true = [200; 50; 0.2];
    end
    if i == round(0.5 * length(t))
        tau_true = [50; -40; -0.5];
    end
    
    % Compute rotation matrix based on current true yaw:
    R = [ cos(X_true(3)), -sin(X_true(3)), 0;
          sin(X_true(3)),  cos(X_true(3)), 0;
          0, 0, 1 ];
      
    % Compute hydrodynamic matrices (based on true state)
    C_M = [ 0, 0, -M(2,2)*X_true(5) - 0.5*(M(2,3)+M(3,2))*X_true(6);
            0, 0, M(1,1)*X_true(4);
            M(2,2)*X_true(5) + 0.5*(M(2,3)+M(3,2))*X_true(6), -M(1,1)*X_true(4), 0 ];
    D_M = -[ Xu + Xuu*abs(X_true(4)), 0, 0;
             0, Yv + Yvv*abs(X_true(5)), Yr;
             0, Nv, Nr + Nrr*abs(X_true(6)) ];
    A_c = [ zeros(3), R;
            zeros(3), -inv(M)*(C_M+D_M) ];
    
    % Discretize continuous dynamics:
    A_d = eye(size(A_c)) + dt * A_c;
    B_d = dt * B_c;
    
    % Update the true state with process noise:
    X_true = A_d * X_true + B_d * tau_true + QF * dt * randn(size(X_true));
    x_true_vec(:, i) = X_true;
end

%% Simulation and AEKF Estimation
for i = 1:length(t)
    % Get current control input from smooth profile
    current_tau_c = [tau_u(i); tau_v(i); tau_r(i)];
    
    % Schedule fault parameter changes
    if i == round(0.5 * length(t))
        theta = [0; 0.041; 0];
    end
    if i == round(0.75 * length(t))
        theta = [0.087; 0.041; 0.2];
    end
    
    % Update system matrices based on current simulation state (eta, nu):
    R = [ cos(eta(3)), -sin(eta(3)), 0;
          sin(eta(3)),  cos(eta(3)), 0;
          0, 0, 1 ];
    C_M = [ 0, 0, -M(2,2)*nu(2) - 0.5*(M(2,3)+M(3,2))*nu(3);
            0, 0, M(1,1)*nu(1);
            M(2,2)*nu(2) + 0.5*(M(2,3)+M(3,2))*nu(3), -M(1,1)*nu(1), 0 ];
    D_M = -[ Xu + Xuu*abs(nu(1)), 0, 0;
             0, Yv + Yvv*abs(nu(2)), Yr;
             0, Nv, Nr + Nrr*abs(nu(3)) ];
    A_c = [ zeros(3), R;
            zeros(3), -inv(M)*(C_M+D_M) ];
    
    % Discretize the system:
    A_d = eye(size(A_c)) + dt * A_c;
    B_d = dt * B_c;
    
    % Fault injection mapping based on current control:
    Psi_d = -B_d * diag(current_tau_c);
    
    % Update the simulation state with fault injection and noise:
    X = A_d * X + B_d * current_tau_c + Psi_d * theta + QF * dt * randn(size(X));
    eta = X(1:3);
    nu  = X(4:6);
    
    % Measurement (full state with noise):
    y = C * X + RF * dt * randn(size(X));
    
    %% Adaptive Extended Kalman Filter (AEKF) Update
    % Shortcut variables for dynamics linearization:
    sig2 = 0.5 * (m23 + m32);
    sig1 = xhat(5) * m22 + xhat(6) * sig2;
    
    % Compute nudot (3x3 matrix for dynamic states)
    nudot = -inv(M) * [ Xu + 2*Xuu*abs(xhat(4)),         xhat(6)*m22,          sig1 + xhat(6)*sig2;
                        -xhat(6)*m11,                     Yv + 2*Yvv*abs(xhat(5)), Yr - xhat(4)*m11;
                         xhat(5)*m11 - sig1,  Nv + xhat(4)*m11 - xhat(4)*m22,  Nr + 2*Nrr*abs(xhat(6)) - xhat(4)*sig2 ];
    
    % Construct the Jacobian FX (6x6)
    % Top block (3x6) based on kinematics:
    top_block = [ 0, 0, -xhat(5)*cos(xhat(3)) - xhat(4)*sin(xhat(3)),  cos(xhat(3)), -sin(xhat(3)), 0;
                  0, 0,  xhat(4)*cos(xhat(3)) - xhat(5)*sin(xhat(3)),  sin(xhat(3)),  cos(xhat(3)), 0;
                  0, 0, 0, 0, 0, 1 ];
    % Bottom block (3x6) with zeros in the first three columns:
    bottom_block = [ zeros(3,3), nudot ];
    FX = A + dt * [ top_block; bottom_block ];
    
    % Kalman Filter prediction update:
    Pmin = FX * Pplus * FX' + QF;
    Sigma = C * Pmin * C' + RF;
    KF = Pmin * C' / Sigma;
    Pplus = (eye(rank(A)) - KF * C) * Pmin;
    
    % Innovation:
    ytilde = y - C * xhat;
    
    % Adaptive updates for noise covariances:
    QF = a * QF + (1 - a) * (KF * (ytilde * ytilde') * KF');
    RF = a * RF + (1 - a) * (ytilde * ytilde' + C * Pmin * C');
    
    % Fault parameter estimation gain update:
    Upsilon = (eye(rank(A)) - KF * C) * FX * UpsilonPlus + (eye(rank(A)) - KF * C) * Psi_d;
    Omega   = C * FX * UpsilonPlus + C * Psi_d;
    Lambda_mat = inv(lambda * Sigma + Omega * S * Omega');
    Gamma = S * Omega' * Lambda_mat;
    S = (1 / lambda) * S - (1 / lambda) * S * Omega' * Lambda_mat * Omega * S;
    UpsilonPlus = Upsilon;
    
    % Update fault parameter estimate:
    thetahat = thetahat + Gamma * (y - C * xhat);
    
    % Update state estimate:
    xhat = A_d * xhat + B_d * current_tau_c + Psi_d * thetahat + QF * dt * randn(size(X)) ...
           + KF * (y - C * xhat) + Upsilon * Gamma * (y - C * xhat);
    
    % Store data for plotting:
    eta_vec(:, i)      = X(1:3);
    nu_vec(:, i)       = X(4:6);
    tau_c_vec(:, i)    = current_tau_c;
    xhat_vec(:, i)     = xhat;
    theta_vec(:, i)    = theta;
    thetahat_vec(:, i) = clamp(thetahat);
end

%% Enhanced Visualization with Descriptive Titles
% Figure 1: Vessel Trajectory
figure(1)
clf;
hold on
plot(eta_vec(1,:), eta_vec(2,:), 'k', 'LineWidth', 3)
plot(xhat_vec(1,:), xhat_vec(2,:), 'r--', 'LineWidth', 3)
plot(eta_vec(1, round(t_control/dt)+1), eta_vec(2, round(t_control/dt)+1), ...
    'gpentagram', 'MarkerSize', 10, 'LineWidth', 2)
title('Vessel Trajectory: True vs Estimated Position')
xlabel('X Position (m)')
ylabel('Y Position (m)')
legend('True Trajectory', 'Estimated Trajectory', 'Control Changes', 'Location', 'best')
grid on

% Figure 2: Body-Fixed Velocities
figure(2)
clf;
subplot(3,1,1)
plot(t, nu_vec(1,:), 'k', t, xhat_vec(4,:), 'r--', 'lineWidth',3)
title('Surge Velocity (u)')
ylabel('m/s')

subplot(3,1,2)
plot(t, nu_vec(2,:), 'k', t, xhat_vec(5,:), 'r--', 'lineWidth',3)
title('Sway Velocity (v)')
ylabel('m/s')

subplot(3,1,3)
plot(t, nu_vec(3,:), 'k', t, xhat_vec(6,:), 'r--', 'lineWidth',3)
title('Yaw Rate (r)')
ylabel('rad/s')
xlabel('Time (s)')
sgtitle('Body-Fixed Velocities: True vs Estimated')
legend('True', 'Estimated', 'Location', 'best')
grid on

% Figure 3: Fault Parameter Estimation
figure(3)
clf;
subplot(3,1,1)
plot(t, theta_vec(1,:), 'k', t, thetahat_vec(1,:), 'r--', 'lineWidth',3)
title('\theta_1 (Surge Actuator Fault)')
ylabel('Magnitude')

subplot(3,1,2)
plot(t, theta_vec(2,:), 'k', t, thetahat_vec(2,:), 'r--', 'lineWidth',3)
title('\theta_2 (Sway Actuator Fault)')
ylabel('Magnitude')

subplot(3,1,3)
plot(t, theta_vec(3,:), 'k', t, thetahat_vec(3,:), 'r--', 'lineWidth',3)
title('\theta_3 (Yaw Actuator Fault)')
ylabel('Magnitude')
xlabel('Time (s)')
sgtitle('Actuator Fault Parameter Estimation')
legend('True Value', 'Estimated Value', 'Location', 'best')
grid on

% Figure 4: Control Inputs
figure(4)
clf;
subplot(3,1,1)
plot(t, tau_c(1,:), 'k', 'LineWidth', 3) 
title('Surge Control Input (\tau_u)')
ylabel('N')

subplot(3,1,2)
plot(t, tau_c(2,:), 'k', 'LineWidth', 3)
title('Sway Control Input (\tau_v)')
ylabel('N')

subplot(3,1,3)
plot(t, tau_c(3,:), 'k', 'LineWidth', 3)
title('Yaw Control Input (\tau_r)')
ylabel('Nm')
xlabel('Time (s)')
sgtitle('Control Inputs Applied to the Vessel')
grid on

%% Clamping Function
function clamped_arr = clamp(arr)
    clamped_arr = arr;
    for i = 1:length(arr)
        x = arr(i);
        if x > 1
            x = 1;
        elseif x < -1
            x = -1;
        end
        clamped_arr(i) = x;
    end
end
