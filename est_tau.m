function tau=est_tau(y,t,y0)
%% Will calculatet the time constant given an array of values and time 
% 1- values
% 2- time 
% 3- initial value can be different than y(1)


% Define the optimization function
fun = @(x) sum((y -y0 * exp(-t/x(1))).^2);

% Set the optimization options
options = optimoptions('fmincon','Display','off');

% Set the initial guess for the time constant
x0 = [0.1];

% Set the lower and upper bounds for the time constant
lb = [0 ];
ub = [1 ];

% Run the optimization
[x,fval,~,output] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);

% Extract the time constant from the optimized parameters
tau = x;

end