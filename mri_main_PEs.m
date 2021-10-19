%% 2nd level parameter estimates
% main experiment

clear
close all
clc

subjID_list = {'sub-MA01'};

% Peaks
peaksMNI = [
    % interaction att x rep
    22 0 52
    -24 -4 54
    -2 14 48
    8 18 40
    -30 26 2
    14 -68 54
    -16 -70 52
    34 -44 46
    -34 -46 46
    ];

peaksNames = {
    'rSFG'
    'lSFG'
    'lCinG'
    'rCinG'
    'lAIns'
    'rSPL'
    'lSPL'
    'rIPS'
    'lIPS'
    };

% Contrasts (check "ve_stim_1level_paramest" contrasts definition)
% 7 = inc
% 8 = con
% 9 = A inv
% 10 = A val
% 11 = V inv
% 12 = V val
% 13 = A inv con
% 14 = A val con
% 15 = V inv con
% 16 = V val con
% 17 = A inv inc
% 18 = A val inc
% 19 = V inv inc
% 20 = V val inc

consToUse = 1:20;

%% Get peaks of interest in matrix coordinates

tmpHdr = spm_vol(fullfile('E:\Data\MAMSI_MRI\', subjID_list{1}, '\derivatives\ve\anl_ve_stim_smooth8\beta_0001.nii'));

m2v=inv(tmpHdr.mat);

for i=1:size(peaksMNI,1)
    peaksVox(i,1:3)=round(peaksMNI(i,:)*m2v(1:3,1:3) + m2v(1:3,4)');
end

clear tmp*

%% Create group-level mean and SEM images
for iCon = 1:length(consToUse)
    
    for iSubj = 1:length(subjID_list)
        subjGLMFolder = fullfile('E:\Data\MAMSI_MRI\', subjID_list{iSubj}, '\derivatives\ve\anl_ve_stim_smooth8');
        tmpHdr = spm_vol(fullfile(subjGLMFolder,...
            ['con_' num2str(consToUse(iCon),'%04i') '.nii']));
        tmpVol = spm_read_vols(tmpHdr);
        allImg(:,:,:,iSubj) = tmpVol;
        clear tmp*
    end
    
    MeanImg{iCon} = mean(allImg,4);
    SEMImg{iCon} = std(allImg,0,4)/sqrt(length(subjID_list));
    
end

%% Plot

% colors
cols.yelAud1  = [165 0 33]/255;
cols.yelAud2  = [255 153 153]/255;
cols.blueVis1 = [153 204 255]/255;
cols.blueVis2 = [0 51 204]/255;

conCols=[cols.blueVis2;cols.blueVis1;cols.blueVis2;cols.blueVis1;...
    cols.yelAud2;cols.yelAud1;cols.yelAud2;cols.yelAud1];

conLineStyles={
    ':'
    '-'
    ':'
    '-'
    };

positionXY = [0, 0, 150, 200];

for iPeak = 1:size(peaksMNI,1)
    
    figure('color', [1 1 1], 'Position', positionXY);
    
    plot_matrix_mean2=[MeanImg{16}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3)) MeanImg{15}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        MeanImg{20}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3)) MeanImg{19}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        MeanImg{13}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3)) MeanImg{14}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        MeanImg{17}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3)) MeanImg{18}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3))];
    
    for i=1:4
        hold on
        if i==1
            k=-0.08;
        elseif i==2
            k=-0.04;
        elseif i==3
            k=0.04;
        elseif i==4
            k=0.08;
        end
        plot((1:2)+k,plot_matrix_mean2(i,:),...
            'Color','k',...
            'LineStyle',conLineStyles{i},...
            'LineWidth',1.5)
    end

    plot_matrix_mean = [MeanImg{16}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        MeanImg{15}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        MeanImg{20}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        MeanImg{19}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        MeanImg{13}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        MeanImg{14}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        MeanImg{17}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        MeanImg{18}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3))];
    plot_matrix_sem = [SEMImg{16}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        SEMImg{15}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        SEMImg{20}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        SEMImg{19}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        SEMImg{13}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        SEMImg{14}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        SEMImg{17}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3));...
        SEMImg{18}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3))];
    
    a=repmat([1 2],1,4);
    
    for i=1:8
        if i==1 || i==2
            k=-0.08;
        elseif i==3 || i==4
            k=-0.04;
        elseif i==5 || i==6
            k=0.04;
        elseif i==7 || i==8
            k=0.08;
        end
        hold on
        line([a(i)+k a(i)+k],[plot_matrix_mean(i)+plot_matrix_sem(i) ...
            plot_matrix_mean(i)-plot_matrix_sem(i)],...
            'Color',conCols(i,:),'LineWidth',1.5);hold on;
        plot(a(i)+k,plot_matrix_mean(i),'o',...
            'Color',conCols(i,:),...
            'LineWidth',1.5,'MarkerSize',4,...
            'MarkerEdgeColor',conCols(i,:),...
            'MarkerFaceColor',conCols(i,:));hold on;
    end
    
    % adjust x and y ticks
    xl = [0.7 2.3]; xlim(xl);
    yl = [0 40]; ylim(yl);
    % adjust x and y ticks
    set(gca,'FontName', 'Helvetica');
    set(gca,'FontSize', 12);
    set(gca,'XTick', [1 2]);
    set(gca,'XTickLabel', []);
    set(gca,'TickLength', [0.02 0.02]);
    set(gca,'LineWidth',1.2)    
    saveas(gcf,fullfile('E:\Data\MAMSI_MRI\group\main',[peaksNames{iPeak} '_inc_inv']),'emf');    
end
