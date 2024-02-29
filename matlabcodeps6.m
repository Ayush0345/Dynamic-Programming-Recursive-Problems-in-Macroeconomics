% ==================================================================
% Example code of value function iteration
% v(k) = max {log(k^alpha + (1-delta)k - k') + beta v(k')}
%
% very crude!
% ==================================================================


% ==================================================================
% Parameters
alpha = 0.36;
beta = 0.95;
delta = 1;  %!! - guess and verify works
nk = 500;  % number of points on the grid of capital
errtol = 10^-6;
% ==================================================================

% Compute exact solution
% v(k) = A log k + B
s = alpha * beta; % savings rate
A = alpha/(1-alpha * beta);
B = (log(1-s) + beta*A* log(s))/(1-beta);


% Find the steady state
kss = ((1/beta - 1 + delta)/alpha)^(1/(alpha-1));


% construct a grid of capital stocks
kgrid = linspace(0.75 * kss, 1.25*kss, nk);
kgrid = kgrid';

vexact = A * log(kgrid) + B;

% ============================= Part a =============================
% Insert code to compute gexact here

gexact = s*kgrid.^alpha;
disp(gexact)

% ============================= Part b =============================
figure(1) 
% Insert code to plot gexact
plot(gexact)
hold on;
plot(kgrid)
hold off;
title("Policy Function v/s Capital")
xlabel("gexact")
ylabel("Capital")
% check command plot

% ==================================================================
% Start the value function iteration
% ==================================================================

Cmatrix = zeros(nk, nk);
Umatrix = zeros(nk, nk);

% NOTE: This is a very inefficient way of computing these matrices! 
% We will learn more efficient ways of doing this!
for i = 1:nk
    for j = 1:nk
        % consumption if current capital is k(i) and future capital is k(j)
        Cmatrix(i,j) = kgrid(i)^alpha + (1-delta)*kgrid(i) - kgrid(j);        
        % utility if current capital is k(i) and future capital is k(j)
        if (Cmatrix(i,j)>0)
            Umatrix(i,j) = log(Cmatrix(i,j));
        else
            Umatrix(i,j) = -inf;
        end
    end
end



% ============================= Part c =============================
% Insert additional code (and modify existing code if necessary) to
% find gn(k)

tic
vguess = zeros(nk, 1);
vnext = zeros(nk, 1);
dif = 1;
while (dif > errtol)
    temp = Umatrix + beta * ones(nk,1) * vguess';
    [vnext, x] = max(temp, [], 2);
    dif = max(abs(vnext - vguess));
    vguess = vnext;
end
toc

gn = zeros(nk, 1);

tic
vguess = zeros(nk, 1);
vnext = zeros(nk, 1);
dif = 1;
while (dif > errtol)
    temp = Umatrix + beta * ones(nk,1) * vguess';
    [vnext, x] = max(temp, [], 2);
    for i = 1:nk
        gn(i) = kgrid(x(i));
    end
    dif = max(abs(vnext - vguess));
    vguess = vnext;
end
toc

disp(gn(500))

% gn = 0.2027
% ========================== Part d ==========================
% In the same plot (use the command hold on) plot gn. 
figure(2)

plot(gexact)
hold on;
plot(gn(500))
hold on;
plot(kgrid)
hold off;
title("Policy Function v/s n-th iteration of Policy Function v/s Capital")
xlabel("gexact & gn")
ylabel("Capital (kgrid)")

% ========================== Part e ==========================
% Write some code that finds the maximum absolute error between 
% gexact and gn. Use the command disp to display it.

error = max(abs(gn - gexact));
disp(error)

% error term = 1.167380578844801e-04
% ========================== Part f ==========================
% Write some code that finds the maximum percentage error between 
% gexact and gn. Use the command disp to display it.

errorpercent = abs(error/gexact);
disp(errorpercent)

% errorpercent = 0.5760



