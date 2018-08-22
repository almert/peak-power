% This code implements a binary search via Monte Carlo simulations for 
% finding the solution to the integral equation in Theorem 3.
% This is for the maximum mutual information setting.

N = 2:35; % Dimensionalities for which to estimate Rmax
T1 = 10000; % Number of data points for the integral over W1
T = 10000; % Number of chi-square samples for each value of W1
epsilon = 0.0001; % Resolution of the value of R
eps_r = 0.0001; % Maximum tolerated residual for the integral equation

R_dim = zeros(1,length(N));
R_res = zeros(1,length(N));
parfor d=1:length(N)
    
    n = N(d);  
    Rb = [sqrt(n),3*sqrt(n)]; % Upper and lower bounds for Rmax
    
    % Binary search is performed in this loop:
    while true
        
        R = mean(Rb);
        integ = zeros(1,T1);
        W1 = linspace(-7,R+7,T1); % Uniform samples from the effective domain of W1
        for i=1:length(W1)
            
            x = W1(i);
            Ws = chi2rnd(n-1,1,T); % Central chi-square with n-1 degrees of freedom           
            W = sqrt(x^2+Ws); % Norm of W
            % W(W<=1e-30) = 1e-30; % To avoid 0/0, if n=1 is considered
            
            % Below are multiple methods to compute the Bessel ratios. To avoid
            % floating point overflows, Steed's or Lentz's methods are advised.
            % This part also integrates over W2:WN.
            %   integ(i)=mean((x./W).*(besseli(n/2,R*W)./besseli(n/2-1,R*W)));
            %   integ(i)=mean((x./W).*arrayfun(@(x) lentzs(n/2,x),R*W));
            integ(i) = mean((x./W).*arrayfun(@(x) steeds(n/2,x),R*W));
            
        end
        % Integrate over W1 and compute the residual of the equation.
        integral = trapz(W1, integ.*(qfunc(W1-R)-qfunc(W1))/R ); 
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
