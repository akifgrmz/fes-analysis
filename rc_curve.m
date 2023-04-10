function RCVarMVC=rc_curve(mean_force,pw_points)

% yRep1=mean_force(1:length(pw_points));
% yRep2=mean_force(length(pw_points)+1:2*length(pw_points));
% yRep3=mean_force(2*length(pw_points)+1:3*length(pw_points));
yRep=mean(mean_force);

init_value = [1, -9, -0.5];

lower_bound = [0, -1000, -1];
upper_bound = [10000, 0, 0];
RCVarMVC = lsqcurvefit(@gompertz,init_value,pw_points,yRep,lower_bound,upper_bound);


end

