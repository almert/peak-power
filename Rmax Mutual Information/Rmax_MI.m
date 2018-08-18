% This code implements a binary search via Monte Carlo simulations for 
% finding the solution to the integral equation in Theorem 3.
% This is for the maximum mutual information setting.

% N = input('Dimension: ');
% T1 = input('How Many Samples in the Integration: ');
% T = input('How samples of Chi-Square to Generate: ');
% epsilon = input('Tolerance value for Rmax: ');
% eps_r = input('Tolerance value for the residual: ');

N = 1:35; % Dimensionalities for which to estimate Rmax
T = 10000; % Number of chi-square samples for each value of W1
T1 = 10000; % Resolution of the integral over W1
epsilon = 0.0001; % Resolution of the value of R
eps_r = 0.0001; % Maximum tolerated residual for the integral equation

X = linspace(-100,100,T1); % Uniform samples from the domain of X.

R_dim = zeros(1,length(N));
R_res = zeros(1,length(N));
parfor d=1:length(N)
    
    n = N(d);  
    Rb = [sqrt(n),3*sqrt(n)]; % Upper and lower bounds for Rmax
    
    % Binary search is performed in this loop:
    while true
        
        R = mean(Rb);
        integ = zeros(1,T1);
        for i=1:length(X)
            
            x = X(i);
            Ws = chi2rnd(n-1,1,T); % Central ci-square with n-1 degrees of freedom
            
            W = sqrt(x^2+Ws);
            
            % Below are multiple methods to compute the Bessel ratios. To avoid
            % floating point overflows, Steed's or Lentz's methods are advised.
            % This part also integrates over W2:WN.
            %   integ(i)=mean((x./W).*(besseli(n/2, r*W )./besseli(n/2-1,r*W)));
            %   integ(i)=mean((x./W).*arrayfun(@(x) lentzs(n/2,x),r*W));
            integ(i) = mean((x./W).*arrayfun(@(x) steeds(n/2,x),R*W));
            
        end
        % Integrate over W1 and compute the residual of the equation.
        integral = trapz(X, integ.*(qfunc(X-R)-qfunc(X))/R ); 
        residual = integral-0.5;
        
        % Check if the residual and the binary search interval are small
        % enough to end the search.
        if abs(residual)<eps_r && Rb(2)-Rb(1)<epsilon
            break;
        elseif residual>0
            Rb(2) = R;
        else
            Rb(1) = R;
        end
    end
    
    R_dim(d) = R; % Save the values of Rmax
    R_res(d) = residual; % Save the final residuals
    
end
save Results_35.mat R_dim R_res N;
