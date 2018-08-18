% This code implements a binary search via Monte Carlo simulations for 
% finding the solution to the second integral equation in the Discussion
% section of our paper. This is for the maximum MMSE setting.

N=1:35; % Dimensionalities for which to estimate Rmax
T1=50000000; % Number of central chi-square samples
T2=50000000; % Number of non-central chi-square samples
epsilon = 0.0001; % Resolution of the value of R
eps_r = 0.0001; % Maximum tolerated residual for the integral equation

Rmaxmax = zeros(1,length(N));
residuals = zeros(1,length(N));
parfor d=1:length(N)
    
    n=N(d); % Dimensionality of the data
    Rb=[sqrt(n), 3*sqrt(n)]; % Upper and lower bounds for Rmax
    
    % Binary search is performed in this loop:
    while true
        
        R=mean(Rb);
        delta= R^2; % Non-centrality parameter
        
        z1 = chi2rnd(n,1,T1); % Central chi-square
        z2 = ncx2rnd(n,delta,1,T2); % Non-central chi-square
        
        % Below are multiple methods to compute the Bessel ratios. To avoid
        % floating point overflows, Steed's or Lentz's methods are advised.
        %   h1=besseli(n/2,sqrt(z1)*R)./besseli(n/2-1,sqrt(z1)*R);
        %   h2=besseli(n/2,sqrt(z2)*R)./besseli(n/2-1,sqrt(z2)*R);
        %   h1 = arrayfun(@(x) lentzs(n/2,x),sqrt(z1)*R);
        %   h2 = arrayfun(@(x) lentzs(n/2,x),sqrt(z2)*R);
        h1 = arrayfun(@(x) steeds(n/2,x),sqrt(z1)*R);
        h2 = arrayfun(@(x) steeds(n/2,x),sqrt(z2)*R);
        
        % Compute the expectation and the residual of the equation.
        integral = mean(h1.^2)+mean(h2.^2);
        residual = integral-1;
        
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
    
    Rmaxmax(d) = R; % Save the values of Rmax
    residuals(d) = residual; % Save the final residuals
    
end
save Results_35.mat Rmaxmax residuals N;
