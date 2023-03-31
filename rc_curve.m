function RCVarMVC=rc_curve(mean_force,pw_points)

yRep1=mean_force(1:length(pw_points));
yRep2=mean_force(length(pw_points)+1:2*length(pw_points));
yRep3=mean_force(2*length(pw_points)+1:3*length(pw_points));
yRep=mean([yRep1;yRep2;yRep3]);

init_value = [10, -90, -0.05];

lower_bound = [0, -1000, -1];
upper_bound = [10000, 0, 0];
RCVarMVC = lsqcurvefit(@gompertz,init_value,pw_points,yRep,lower_bound,upper_bound);


end

