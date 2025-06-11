%% Data Inject
clear all

TestFolders=["nov21_24" "nov27_24" "dec12_24"];
FileNames=strings(1,length(TestFolders))+"expsave";

for iFile=1:length(TestFolders)
    ExpStruct=sprintf('%s',char(FileNames(iFile)));
    str=sprintf('%s/%s',char(TestFolders(iFile)),ExpStruct);
    M=load (str);
    FldNames = fieldnames(M);

    for iField=1:length(FldNames)
        S.(TestFolders(iFile))=M.(FldNames{iField});
    end
end



%% Sigmoid f (Tracking Error)
% define parameters and functions
% func_h: hand error sigmoid function
% func_e: tracking error function
close all
theta_np = 1;  % Example coefficient
E = 1;         % Example energy value
thres_c = 0.1;         % Threshold for e_tr
thres_h = 0.4;   
thres_E = 0.5;

% Threshold for Error_hand
kf_values=[1 3  5  10 25]; 
kh_values=[3 5  10 20 40] ; 
kE_values=[5 10 20 50 100];
% e_tr_values = linspace(-1, 1, 100);
% e_h_values = linspace(0, 1, 100);

% kf=[  0.5 1.5  2  3 5 10]; 
% kh=[  1   1.2  2  3 5]; 
w_h=0.5;
w_f=0.5;
% thres_c=0.3;
% h=0.5;
% % calculate the sigmoids
e_tr_values=linspace(-0.5, 1.2, 200);
e_h_values=linspace(-0.8, 0.8, 200);
Ef_values=linspace(0, 1, 100);
sigmoid_f=zeros(length(e_tr_values),length(kf_values));
sigmoid_h=zeros(length(e_h_values),length(kh_values));
sigmoid_E=zeros(length(Ef_values),length(kE_values));
function f_e = func_f( e_tr, c, k1)
    % Smooth transition function for e_tr
    % e_tr=e_tr+c;
    % f_e = (1 / (1 + exp(-k1 * (e_tr - c)))) * min(1,(e_tr / c)^3);
    f_e = (1 / (1 + exp(-k1 * (e_tr - c))));
    
    f_e(f_e<0)=0;
    % f_e=f_e+c;
end

function f_h = func_h( error, h, k)
    % error=error-h;
    error=abs(error);
    f_h = (1 / (1 + exp(k * (error-h ))));
end

function f_E = func_E( E, ef, ke)
    % error=error-h;
    E=abs(E);
    f_E = (1 / (1 + exp(ke * (ef-E ))));
end


for j=1:length(kf_values)
    for i = 1:length(e_tr_values)
        sigmoid_f(i,j) = func_f(e_tr_values(i),thres_c,kf_values(j));
    end
end

for j=1:length(kh_values)
    for i = 1:length(e_h_values)
        sigmoid_h(i,j) = func_h(e_h_values(i),thres_h,kh_values(j));
    end
end

for j=1:length(kE_values)
    for i = 1:length(Ef_values)
        sigmoid_E(i,j) = func_E(Ef_values(i),thres_E,kE_values(j));
    end
end

% sigmoid_f(sigmoid_f < 0) = 0;        
figure(10)
subplot(3,1,1)
hold on

legend_tr = cell(1, length(kf_values));
for j = 1:length(kf_values)
    plot(e_tr_values, sigmoid_f(:,j),'LineWidth',2)
    legend_tr{j} = sprintf('k_f = %.1f', kf_values(j));
end
plot([ 0 0],[-1 1],'--k','LineWidth',0.01)
% plot([ thres_c thres_c],[0 1],'--k','LineWidth',0.01)

grid on 
xlabel('e_{tr}')
ylabel('sigmoid f')
ylim([-0.2 1])
xlim([-0.2 0.5])
legend (legend_tr)

subplot(3,1,2)
hold on
legend_h = cell(1, length(kh_values));
for j = 1:length(kh_values)
    plot(e_h_values, sigmoid_h(:,j),'LineWidth',2)
    legend_h{j} = sprintf('k_h = %.1f', kh_values(j));
end
grid on 
xlabel('e_{hand}')
ylabel('sigmoid h')
ylim([0 1])
legend (legend_h)

subplot(3,1,3)
hold on
legend_E = cell(1, length(kE_values));
for j = 1:length(kE_values)
    plot(Ef_values, sigmoid_E(:,j),'LineWidth',2)
    legend_E{j} = sprintf('k_h = %.1f', kE_values(j));
end
grid on 
xlabel('Effort')
ylabel('sigmoid E')
ylim([0 1])
legend (legend_E)
%%
% theta_np = 1;  % Example coefficient
% E = 1;         % Example energy value
% thres_c = 0.1;         % Threshold for e_tr
% thres_h = 0.4;         % Threshold for Error_hand
% kf_values=[  0.5 1.5  2  3 5 10]; 
% kh_values=[  1   1.2  2  3 5]; 
% e_tr_values = linspace(-1, 1, 100);
% e_h_values = linspace(0, 1, 100);

S_values = zeros(length(e_tr_values), length(e_h_values));
S_values2 = zeros(length(e_tr_values), length(e_h_values));

num_kf=length(kf_values);
num_kh=length(kh_values);
num_kE=length(kE_values);
plot_idx=1;
S_idx=1;
S_mat=zeros(num_kf*num_kh*length(e_tr_values)*length(e_h_values),6);
for i_kf=1:num_kf
    kf=kf_values(i_kf);

    for i_kh=1:num_kh
        kh=kh_values(i_kh);
        for i_E=1:num_kE
            kE=kE_values(i_E); 
            for i_etr = 1:length(e_tr_values)
                for i_hand = 1:length(e_h_values)
                    sigmoid_f(i_etr,i_hand) = func_f(e_tr_values(i_etr),thres_c,kf);
                    sigmoid_h(i_etr,i_hand) = func_h(e_h_values(i_hand),thres_h,kh);  

                    for i_effort = 1:length(Ef_values)
                        sigmoid_E(i_etr,i_hand,i_effort) = func_E(Ef_values(i_effort),thres_E,kE);
    
                        S_values(i_etr,i_hand,i_effort) =sigmoid_E(i_etr,i_hand,i_effort)*theta_np*sigmoid_f(i_etr,i_hand) *sigmoid_h(i_etr,i_hand) ;
                        S_values2(i_etr,i_hand,i_effort) =sigmoid_E(i_etr,i_hand,i_effort)*theta_np*(w_f*sigmoid_f(i_etr,i_hand) + w_h*sigmoid_h(i_etr,i_hand)) ;
                        % S_val_E(i_etr,i_effort) =Ef_values(i_effort)*theta_np*(w_f*sigmoid_f(i_etr,i_hand) + w_h*1) ;
                        
                        % S_mat(S_idx,:)=[S_values(i_etr,i_hand,i_effort) S_values2(i_etr,i_hand,i_effort) kf kh  thres_c thres_h] ;
                        S_idx=S_idx+1;
                    end
                end
            end
                                            figure(13);
            set(gcf,'Position', [100, 100, 1400, 1000]);
            subplot(num_kf,num_kh,plot_idx)
            imagesc(e_h_values, e_tr_values, S_values(:,:,end));
            title(sprintf("k_f=%.2f, k_h=%.2f",kf,kh))
            xlabel('e_{Hand}')
            ylabel('e_{track}')
            caxis([0 1]); % Fix color axis
            % 
                    figure(14);
            set(gcf,'Position', [100, 100, 1400, 1000]);
            subplot(num_kf,num_kh,plot_idx)
            imagesc(e_h_values, e_tr_values, S_values2(:,:,end));
            title(sprintf("k_f=%.2f, k_h=%.2f",kf,kh))
            xlabel('e_{Hand}')
            ylabel('e_{track}')
            caxis([0 1]); % Fix color axis
    
            figure(15);
            set(gcf,'Position', [100, 100, 1400, 1000]);
            subplot(num_kE,num_kf,plot_idx)
            imagesc(Ef_values,e_tr_values,squeeze(S_values2(:,100,:)));
            title(sprintf("k_f=%.2f, k_E=%.2f",kf,kE))
            xlabel('e_{Effort}')
            ylabel('e_{track}')
            caxis([0 1]); % Fix color axis
            plot_idx=plot_idx+1;
        end
    end
    colorbar
end

%% smooth cont plotting
close all
Trials=[12]; %% 12 and 15 is the good examples 
for iTest=1:length(TestFolders)
    TestLabel=TestFolders(iTest);
    ExpLabel=S.(TestLabel).ExpLabels(6);
    EffortFESTable=S.(TestLabel).(ExpLabel).EffortFESTable;
    CompTrialNum=sum(EffortFESTable(:,end));
    TrialNum=length(EffortFESTable(:,end));

    for iTrial=1:length(Trials)
        TrialLabel=sprintf("Trial_%d",Trials(iTrial));
        
         % rows are frames, columns are samples in a frame
        ign=9;  
        FiltvEMG_f=S.(TestLabel).(ExpLabel).(TrialLabel).FiltvEMG_f(ign:end,:);
        FiltmWaves_f=S.(TestLabel).(ExpLabel).(TrialLabel).FiltmWaves_f(ign:end,:);
        Trig_f=S.(TestLabel).(ExpLabel).(TrialLabel).Trig_f(ign:end,:);
        BlankedEMG_f=S.(TestLabel).(ExpLabel).(TrialLabel).BlankedEMG_f(ign:end,:);
        Dropped_f=S.(TestLabel).(ExpLabel).(TrialLabel).DroppedEMG_f(ign:end,:);
        DroppedTrig_f=S.(TestLabel).(ExpLabel).(TrialLabel).DrpdTrig_f(ign:end,:);
        PW_f=S.(TestLabel).(ExpLabel).(TrialLabel).PW_f(ign:end,:);
        
        Time_f=S.(TestLabel).(ExpLabel).(TrialLabel).Time_f(ign:end,:);                
        Tm_f=S.(TestLabel).(ExpLabel).(TrialLabel).Tm_f(ign:end,:);
        Time_f=Tm_f;
        [r,c]=size(Time_f);
        
        FiltvEMG_f=reshape(FiltvEMG_f',1,[]);
        Time_f=reshape(Time_f',1,[]);
        FiltmWaves_f=reshape(FiltmWaves_f',1,[]);
        Trig_f=reshape(Trig_f',1,[]);
        BlankedEMG_f=reshape(BlankedEMG_f',1,[]);                
        Dropped_f=reshape(Dropped_f',1,[]);
        DroppedTrig_f=reshape(DroppedTrig_f',1,[]);
        PW_f=reshape(PW_f',1,[]);

        Hands=S.(TestLabel).(ExpLabel).(TrialLabel).data(:,4:5);
        PW=S.(TestLabel).(ExpLabel).(TrialLabel).data(:,6);
        Target=S.(TestLabel).(ExpLabel).(TrialLabel).data(:,21);
        Stim=S.(TestLabel).(ExpLabel).(TrialLabel).data(:,22);
        Adap_filt=S.(TestLabel).(ExpLabel).(TrialLabel).data(:,23);
        Time=S.(TestLabel).(ExpLabel).(TrialLabel).time;

                %%Smooth cont
        Error_t=Target-Hands(:,1) +0.1;
        Error_h=abs(Hands(:,2)-Hands(:,1));
        Effort=Adap_filt;
        theta_np=Hands(:,2);

        % Define parameters
        c = 0.2;         % Threshold for e_tr
        h = 0.4;         % Threshold for Error_hand
        k1 = 1;        % Smoothness factor for e_tr
        k2 = 1;        % Smoothness factor for Error_hand

        % Store results
        S_values = zeros(length(Error_t),1);

        % Compute S_i for different values
        for i = 1:length(Error_t)
                S_values(i) = 3.5*smooth_S(theta_np(i), Effort(i), Error_t(i), c, Error_h(i), h, k1, k2)+0.08333333;
        end
                S_values(S_values < 0.08333333) = 0.08333333;        
        % S values Map
        % Generate sample values for e_tr and Error_hand
        e_tr_values = linspace(-1, 1, 100);
        Error_hand_values = linspace(0, 2, 100);
        the_np=0.3;
        E=0.5;
        % Store results
        S_val = zeros(length(e_tr_values), length(Error_hand_values));
        
        % Compute S_i for different values
        for i = 1:length(e_tr_values)
            for j = 1:length(Error_hand_values)
                S_val(i, j) = smooth_S(the_np, E, e_tr_values(i), c, Error_hand_values(j), h, k1, k2);
            end
        end
        
        % Plot the results
        figure(5);
        imagesc(Error_hand_values, e_tr_values, S_val);
        colorbar;
        xlabel('Error_{hand}');
        ylabel('e_{tr}');
        title('Smooth Transition of S_i');


        ind=100:r*c-100;
        f=figure(iTest);
        % f.Position = [800 300 800 1500 ];    
        % subplot(length(Trials),2,2*iTrial)

        subplot(2,1,1)        
        plot(Time,Hands(:,1),'LineWidth',1,"DisplayName","Paretic")
        hold on
        plot(Time,Hands(:,2),'LineWidth',1,"DisplayName","NonParetic")

        % plot(Time,MAVs(:,1),'LineWidth',1,"DisplayName","vEMG MAV")
        % plot(Time,MAVs(:,2),'LineWidth',1,"DisplayName","mWave MAV")
        % plot(Time,MAVs(:,3),'LineWidth',1,"DisplayName","Unfiltered")
        % plot(Time,MAVs(:,4),'LineWidth',1,"DisplayName","Dropped")                
        plot(Time,Target,'--','LineWidth',2,"DisplayName","Target")
        plot(Time,S_values,'LineWidth',1,"DisplayName","Updated Cont.")
        plot(Time,Adap_filt,'LineWidth',1,"DisplayName","Effort Est.")
        plot(Time,PW,'LineWidth',1,"DisplayName","Previous Cont.")        
        ylim([-.1 1.1])        
        legend
        xlabel('Time (s)')     
        ylabel('Norm. Tracking Signals')            
        

%         % subplot(length(Trials),2,2*iTrial-1)
%         subplot(2,1,2)        
%         % plot(Time_f(ind),Trig_f(ind)/200    ,'LineWidth',0.5,"DisplayName","Trig") 
%         hold on        
%         plot(Time_f(ind),BlankedEMG_f(ind),'LineWidth',1,"DisplayName","BlankedEMG")  
%         plot(Time_f(ind),FiltmWaves_f(ind),'LineWidth',1,"DisplayName","FiltmWaves")        
%         plot(Time_f(ind),FiltvEMG_f(ind),'LineWidth',1,"DisplayName","FiltvEMG")
%         % plot(Time_f(ind),DroppedTrig_f(ind)/100,'LineWidth',0.5,"DisplayName","DroppedTrig")
%         % plot(Time_f(ind),Dropped_f(ind),'LineWidth',1,"DisplayName","DroppedEMG")
% %         plot(Time_f(ind),PW_f(ind),'LineWidth',1,"DisplayName","PW")
% %         legend({'Trigger(a.u.)','Raw EMG (a.u.)'})
%         legend
%         title(TrialLabel)
%         xlabel('Time (s)')
%         ylim([-.1 .1])
% 
    
        
    end
end



%% smooth controller (NO USE)

% Define parameters
theta_np = 1;  % Example coefficient
E = 5;         % Example energy value
c = 0.1;         % Threshold for e_tr
h = 0.4;         % Threshold for Error_hand
k1 = 2;        % Smoothness factor for e_tr
k2 = 0.5;        % Smoothness factor for Error_hand
% Generate sample values for e_tr and Error_hand
e_tr_values = linspace(-1, 1, 100);
Error_hand_values = linspace(0, 1, 100);

% Store results
S_values = zeros(length(e_tr_values), length(Error_hand_values));

% Compute S_i for different values
for i = 1:length(e_tr_values)
    for j = 1:length(Error_hand_values)
        S_values(i, j) = smooth_S(theta_np, E, e_tr_values(i), c, Error_hand_values(j), h, k1, k2);
    end
end

% Plot the results
figure;
imagesc(Error_hand_values, e_tr_values, S_values);
colorbar;
xlabel('Error_{hand}');
ylabel('e_{tr}');
title('Smooth Transition of S_i');

% Define the function smooth_S
function S_i = smooth_S(theta_np, E, e_tr, c, Error_hand, h, k1, k2)
    % Smooth transition function for e_tr
    f_e = (1 / (1 + exp(-k1 * (e_tr - c)))) * min(1, e_tr / c);
    % f_e = (1 / (1 + exp(-k1 * e_tr))) * (1 - exp(-c * e_tr));
    
    % Ensure f_e is non-negative
    % f_e = max(f_e, 0);

    % Smooth transition function for Error_hand
    g_E = 1 / (1 + exp(-k2 * (Error_hand - h)));
    
    % Compute smooth S_i
    S_i = theta_np * E * f_e * g_E;
end

