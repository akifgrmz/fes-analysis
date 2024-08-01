
close all
k=1;
V=[ 0.2 0.4 0.6 0.8];
S=0.01:0.01:1;
O_v=1.5;
for i=1:length(V)
    T=V(i)+S-k.*V(i).*S;

    V_o=V(i).*T./(V(i)+S)*O_v;
    S_o=S.*T./(V(i)+S);
    
    figure(1)
    subplot(2,2,1)
    plot(S,V_o)
    hold on
    grid on
    subplot(2,2,3)
    plot(S,S_o)
    hold on
    grid on

end


k_s=1.5;
k_v=.6;

V=[ 0.2 0.4 0.6 0.8];
S=0.01:0.01:1;
for i=1:length(V)

    T_s=V(i)+S-k_s.*V(i).*S;
    T_v=V(i)+S-k_v.*V(i).*S;

    V_o=V(i).*T_v./(V(i)+S);
    S_o=S.*T_s./(V(i)+S);
    
    figure(1)
    subplot(2,2,2)
    plot(S,V_o,'LineWidth',2,'DisplayName',sprintf('V=%d',V(i)))
    hold on
    grid on
end

S=[ 0.2 0.4 0.6 0.8];
V=0.01:0.01:1;
for i=1:length(S)


    T_s=V+S(i)-k_s.*V.*S(i);
    T_v=V+S(i)-k_v.*V.*S(i);

    V_o=V.*T_v./(V+S(i));
    S_o=S(i).*T_s./(V+S(i));
    subplot(2,2,4)
    plot(V,S_o,'LineWidth',2,'DisplayName',sprintf('V=%d',V(i)))
    hold on
    grid on

end