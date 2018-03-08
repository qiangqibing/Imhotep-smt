%% Generate the system for Test Case 1
clc; close all; clear all;
rand('state',0);
randn('state',0);

for n = [10, 25, 50, 75, 100, 150]
    clear X  Y  E O Y_ETPG Ablkdiag Cblkdiag
    p           = 20;
    tau         = n;
    % Generate a system with a random A matrix
    A           = full(sprand(n,n,0.3)); %randn(n,n);
    % Make sure A has spectral radius 1 (otherwise A^k will be either very
    % large or very small for k large)
    A           = A/(max(abs(eig(A))) + 0.1);
    % The 'C' matrix of the system
    C           = full(sprand(p,n,0.2)); %randn(p,n);
    Ts          = 0.1;
    
    sys         = ss(A,zeros(n,1),C,0, Ts);
    x0          = randn(n,1);
    
    attackpower = 20;   % Magnitude of the attacks (i.e., norm of the attack vector)
    max_s       = floor(p/2-1)-1-3; % limit to only 5 sensors to be attacked
    s           = max_s;
    
    % Choose a random attacking set K of size qs
    per         = randperm(p);
    K           = per(1:s);
    
    % Choose an initial condition
    x = x0;
    Y = [];
    
    for t = 1 : 1 : tau
        t
        % Generate a random attack vector supported on K
        a       = zeros(p,1);
        a(K)    = attackpower*randn(length(K),1);
        % The measurement is y=C*x+a
        y       = C*x + a;
        % Update the arrays X,Y,E
        Y       = [Y y];
        
        x       = A*x;
    end
    
    save(['./Test1_states/test_n' num2str(n) '_p' num2str(p)]);
end


%% Housekeeping
clc; close all; clear all;  rand('state',0);    randn('state',0);
%imhotepSMTPath = '../../';  %path to Imhotep-SMT
%addpath(imhotepSMTPath);


%% In this test, the number of sensors is fixed to 20 and we increase the
% number of states from 10 - 150. The systems are pre-simulated and the
% inputs and outputs are recorded.

% Generate the system
format longg
test_counter        = 0;
TimeSpent_SMT_test1 = [];
TimeSpent_CVX_test1 = [];

RelativeEstimationError_SMT_test1 = [];
RelativeEstimationError_CVX_test1 = [];
p = 20;     % the number of sensors is fixed to 20, 5 of them under attack
for n = [10, 25, 50, 75, 100, 150]
    % load the system
    load(['./Test1_states/test_n' num2str(n) '_p' num2str(p)]);
    test_counter = test_counter + 1;
    
    disp(['==== Running Test number 1 with n = ' num2str(n) ' and p = ' num2str(p) ' ====']);
    disp(' '); disp(' ');
    % configure ImhotepSMT
    noiseBound                  = zeros(p,1);       % noise bounds
    max_sensors_under_attack    = int8(8);          % maximum sensors under attack
    
    safe_sensors                = [];
    
    smt = ImhotepSMT();
    smt.init(sys, max_sensors_under_attack, safe_sensors, noiseBound );
    
    estimatedAttackedSensorIndex={};
    for t = 1 : n
        % load the outputs from collected data
        y = Y(:,t);
        tic;
        [xhat, b] = smt.addInputsOutputs(0, y);
        time = toc;
        estimatedAttackedSensorIndex{t}=b;
       % disp(['Execution time = ' num2str(time) ' sec']);
        TimeSpent_SMT_test1(test_counter) = time;
        RelativeEstimationError_SMT_test1(test_counter)=norm(xhat-x0)/norm(x0);
    end
    disp(' ');
    
    
    % below for l1 convex optimization
    
    internal_n                       = size(sys.A,1);
    internal_p                       = size(sys.C,1);
    internal_m                       = size(sys.B,2);
    internal_s                       = double(max_sensors_under_attack);
    internal_tau                     = internal_n;
    
    for sensorIndex = 1 : internal_p
        internal_O{sensorIndex}          = obsv(sys.A, sys.C(sensorIndex,:));
        internal_Y{sensorIndex}          = zeros(tau,1);
    end
    
    for sensorIndex = 1 : internal_p
        internal_Y{sensorIndex}     = Y(sensorIndex, :)';
    end
    
    stacked_Y=[];
    stacked_O=[];
    for sensorIndex = 1 : internal_p
        
        stacked_Y=[stacked_Y , internal_Y{sensorIndex}];
        stacked_O=[stacked_O , internal_O{sensorIndex}];
        
    end
    
    % start: below is for ECOS Julia
    
    catenate_Y=zeros(internal_n, internal_p);
    catenate_O=zeros(internal_n, internal_n, internal_p);
    
    for sensorIndex=1: internal_p
        catenate_Y(:,sensorIndex)=internal_Y{sensorIndex};
        catenate_O(:,:,sensorIndex)=internal_O{sensorIndex};
    end
    save(['./Test1_states/test_n' num2str(n) '_p' num2str(p)]);
    save(['./dataForECOSJulia/ECOSJulia_n' num2str(internal_n) '_p' num2str(internal_p)], 'internal_n', 'internal_p', 'catenate_Y', 'catenate_O')
    
    % end
    
    tic
    cvx_begin
    variable state(internal_n)
    variable attack(internal_n, internal_p)
    
    minimize sum(norms(attack, 2))
    
    subject to
    
    for sensorIndex= 1: internal_p
        internal_Y{sensorIndex} == internal_O{sensorIndex}*state+attack(:,sensorIndex)
    end
    
    cvx_end
    
    time= toc;
    TimeSpent_CVX_test1(test_counter) = time;
    RelativeEstimationError_CVX_test1(test_counter)=norm(state-x0)/norm(x0);
    
end

% Plot figure
figure;
plot([10, 25, 50, 75, 100, 150], TimeSpent_SMT_test1,'LineWidth',3)
hold on

plot([10, 25, 50, 75, 100, 150], TimeSpent_CVX_test1,'LineWidth',3)
set(gca,'FontSize',30);
xlabel('number of states');
ylabel('Time (sec)');



figure;
plot([10, 25, 50, 75, 100, 150], RelativeEstimationError_SMT_test1,'LineWidth',3)
hold on

plot([10, 25, 50, 75, 100, 150], RelativeEstimationError_CVX_test1,'LineWidth',3)
set(gca,'FontSize',30);
xlabel('number of states');
ylabel('Relative Estimation Error (sec)');



