function [shearflows_numeric] = shearanalysis(A, Ycent, Zcent, Fo, Iyy)

% Getting Number of Lumps 
n = numel(A);

% Creating the setup matrix for the shearflows in terms of q1 and distances
syms q1 real
shearflows_symbolic = sym(zeros(1,n));
distances = sym(zeros(1,n));
shearflows_symbolic(1) = q1; % q1 is know from the first point


% Loop to find all shear flows in terms of q1 using a correction factor,
% which is the change in the shear flow from lump to lump. 
for j = 2:n
    correction = 0;
    for k = 1:(j-1)
        correction = correction + Zcent(k) * A(k+1);
    end
    shearflows_symbolic(j) = q1 - (Fo/Iyy) * correction;
end



% Step 2: Calculate distances between points (s) for the moment equalibrium

for j = 1:n-1
    distances(j) = sqrt( (Ycent(j) - Ycent(j+1))^2 + (Zcent(j) - Zcent(j+1))^2 );
end
distances(n) = sqrt( (Ycent(n) - Ycent(1))^2 + (Zcent(n) - Zcent(1))^2 );



%Moment Arm from point 1 to q1 to qn-1
for j = 1:n-1
   
    Yavg = (Ycent(j) + Ycent(j+1)) / 2; 
    Zavg = (Zcent(j) + Zcent(j+1)) / 2;
    
    moment_arms(j) = (Zavg - Zcent(1)) * (Ycent(j+1) - Ycent(j)) - (Yavg - Ycent(1)) * (Zcent(j+1) - Zcent(j));
end


% Amoment Arm from point 1 to from qn to q1
Yavg = (Ycent(n) + Ycent(1)) / 2;
Zavg = (Zcent(n) + Zcent(1)) / 2;
moment_arms(n) = (Zavg - Zcent(1)) * (Ycent(1) - Ycent(n)) - (Yavg - Ycent(1)) * (Zcent(1) - Zcent(n));

% Solving the equation of moment equalibrium 
moment_eqn = 0;
for j = 1:n
    moment_eqn = moment_eqn + shearflows_symbolic(j) * moment_arms(j);
end

% Solving moment equalibrium equation for q1
q1_value = solve(moment_eqn == 0, q1);

% Using q1 to calculate all shear flows 
shearflows_numeric = double(subs(shearflows_symbolic, q1, q1_value));


end