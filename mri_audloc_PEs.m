%% 2nd level parameter estimates
% auditory localiser

clear
close all
clc

subjID_list = {'sub-MA01'};

% Peaks
peaksMNI = [
    18 -52 -24 % rCer (t>b)
    -50 -30 10 % lPT (t>b)
    64 -34 16 % rPT (t>b)
    64 -26 22 % rPOp (t>b)
    -36 -28 50 % lPostG (t>b)
    -36 -12 64 % lPreG (t>b)
    26 -16 46 % rSFG (t>b)
    -22 -16 52 % lSFG (t>b)
    
    26 -16 74 % rPreG (L>R)
    34 -22 6 % rHG (L>R)
    52 -24 22 % rPOp (L>R) ***
    50 -32 18 % rPT (L>R) ***
    52 -22 4 % rPT (L>R)
    22 -52 72 % rSPL (L>R) ***
    
    52 20 28 % rIFG (R>L)
    -36 -28 52 % lPostG (R>L)
    -50 -32 8 % lPT (R>L) ***
    -60 -36 14 % lPOp (R>L) ***
    -8 -74 38 % lPreCu (R>L)
    -38 -54 58 % lSPL (R>L) ***
    38 -68 48 % rSPL (R>L)
    ];

peaksNames = {
    'rCer(t-b)'
    'lPT(t-b)'
    'rPT(t-b)'
    'rPOp(t-b)'
    'lPostG(t-b)'
    'lPreG(t-b)'
    'rSFG(t-b)'
    'lSFG(t-b)'
    
    'rPreG(L-R)'
    'rHG(L-R)'
    'rPOp(L-R)'
    'rPT(L-R)'
    'rPT(L-R)'
    'rSPL(L-R)'
    
    'rIFG(R-L)'
    'lPostG(R-L)'
    'lPT(R-L)'
    'lPOp(R-L)'
    'lPreCu(R-L)'
    'lSPL(R-L)'
    'rSPL(R-L)'
    };

% Contrasts (check "audloc_stim_1level_paramest" contrasts definition)
consToUse = [1 2 3];

%% Get peaks of interest in matrix coordinates

tmpHdr = spm_vol(fullfile('E:\Data\MAMSI_MRI\', subjID_list{1}, '\derivatives\audloc\anl_audloc_smooth8\beta_0001.nii'));

m2v=inv(tmpHdr.mat);

for i=1:size(peaksMNI,1)
    peaksVox(i,1:3)=round(peaksMNI(i,:)*m2v(1:3,1:3) + m2v(1:3,4)');
end

clear tmp*

%% Create group-level mean and SEM images
for iCon = 1:length(consToUse)
    
    for iSubj = 1:length(subjID_list)
        subjGLMFolder = fullfile('E:\Data\MAMSI_MRI\', subjID_list{iSubj}, '\derivatives\audloc\anl_audloc_smooth8');
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

conCols = [
    0.2 0.2 0.2
    0.5 0.5 0.5
    0.8 0.8 0.8
    ];

conLineStyles = {
'-'
'-'
'-'
};

conNames = {
    'L'
    'C'
    'R'    
    };

positionXY = [0, 0, 80, 150];

for iPeak = 1:size(peaksMNI,1)
    %     for iPeak = 1
    
    figure('color', [1 1 1], 'Position', positionXY);
    hold on
    
    for iCon = 1:length(conNames)
        bar(iCon,MeanImg{iCon}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3)),...
            0.7,'FaceColor',conCols(iCon,:),'LineStyle','none')
        line([iCon iCon],[...
            MeanImg{iCon}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3))+...
            SEMImg{iCon}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3)) ...
            MeanImg{iCon}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3))-...
            SEMImg{iCon}(peaksVox(iPeak,1),peaksVox(iPeak,2),peaksVox(iPeak,3))],...
            'Color',[0 0 0],'LineWidth',1.2)
    end
    
    title([peaksNames{iPeak} ' ' num2str(peaksMNI(iPeak,:))])
    
    set(gca,'FontName', 'Helvetica');
    set(gca,'FontSize', 6);
    set(gca,'XTick',[]);
    set(gca,'TickLength', [0.03 0.03]);
    xlim([0.2 3.8]);
    set(gca, 'LineWidth', 1.2)
    ylim([-4 4]);
    set(gca, 'YTick', (-4:2:4));
    set(gca,'XColor','none')
    saveas(gcf,fullfile('E:\Data\MAMSI_MRI\group\audloc',peaksNames{iPeak}),'emf');
end
