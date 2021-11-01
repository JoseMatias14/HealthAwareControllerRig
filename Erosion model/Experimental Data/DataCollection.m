% This code computes the change in the probe area along the experiment and compares it with the pressure drop. 


% Author: Jose Otavio Matias
% Work address
% email: jose.o.a.matias@ntnu.no
% May 2020; Last revision: 2020-05-26
clear 
close all
clc
% 


%% Parameters tuning 
%number of wells being used in the experiment
par.nw = 3;

% Sampling times 
par.dtMeasu = 1; %[s]

% Sampling times 
par.dtImage = 15; %[s]

% Sampling times 
par.Image0 = 13; %[s] experiments starts 11:30:03 | Image starts 11:30:15

%threshold for the color when identifying the band in the cropped image
par.colorThreshold =  110;

%filter level 
par.fL = 10;

% for saving
filename = 'Experiment_2';

% File with the data
name = [pwd,'/2021-07-08_112958_Exp2.txt'];

%%  getting the images
par.currentdirectory = pwd;
par.nfolder = dir([par.currentdirectory,'/Images/well_1/*.png']);%pictures taken with png format
par.ni = size(par.nfolder,1); %number of images (all cameras have the same number of pictures

% arrays for storing the probe area
area = zeros(par.nw,par.ni);


%% reading data 
opts = delimitedTextImportOptions("NumVariables", 36);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Date", "Time", "Relativetimes", "FIC104setptslmin", "FIC104slmin", "FIC105setptslmin", "FIC105slmin", "FIC106setptslmin", "FIC106slmin", "CV101setptlmin", "FI101lmin", "CV102setptlmin", "FI102lmin", "CV103setptlmin", "FI103lmin", "dP101mbarD", "dP102mbarD", "dP103mbarD", "CV107setptmbarG", "PI101mbarG", "CV108setptmbarG", "PI102mbarG", "CV109setptmbarG", "PI103mbarG", "PumpoutputpressuresetptbarG", "PI104barG", "TI101C", "TI102C", "TI103C", "CV101currentA", "CV102currentA", "CV103currentA", "CV107currentA", "CV108currentA", "CV109currentA", "PumpctrlcurrentA"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Date", "Time"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["Date", "Time", "Relativetimes", "FIC104setptslmin", "FIC104slmin", "FIC105setptslmin", "FIC105slmin", "FIC106setptslmin", "FIC106slmin", "CV101setptlmin", "FI101lmin", "CV102setptlmin", "FI102lmin", "CV103setptlmin", "FI103lmin", "dP101mbarD", "dP102mbarD", "dP103mbarD", "CV107setptmbarG", "PI101mbarG", "CV108setptmbarG", "PI102mbarG", "CV109setptmbarG", "PI103mbarG", "PumpoutputpressuresetptbarG", "PI104barG", "TI101C", "TI102C", "TI103C", "CV101currentA", "CV102currentA", "CV103currentA", "CV107currentA", "CV108currentA", "CV109currentA", "PumpctrlcurrentA"], "DecimalSeparator", ",");
opts = setvaropts(opts, ["Date", "Time"], "ThousandsSeparator", ".");

expData = table2array(readtable(name, opts));
expTime = size(expData,1);

for cc = 1:par.nw
    data{cc}.flowrate = expData(:, 11 + 2*(cc - 1));
    data{cc}.deltaP = expData(:, 16 + (cc - 1));
    data{cc}.temperature = expData(:, 27 + (cc - 1));
    data{cc}.ptop = expData(:, 20 + 2*(cc - 1));
    data{cc}.ppump = expData(:, 26);
    data{cc}.pumpRotation = (12 + (92 - 12)*(expData(:, 36) - 0.004)./(0.02 - 0.004))./100; % converting from mA to [0 - 1]
    data{cc}.valveOpen = (expData(:, 30 + (cc - 1)) - 0.004)./(0.02 - 0.004); % converting from mA to [0 - 1]
    
    data{cc}.units = {'Q [L/min] ','dP [mbar]','T [oC]','P_{top} [mbar g]','P_{pump} [bar g]','v_{open} [0-1]','P_{rot} [0-1]'}; 
end


%% Computing area
% loop to compute the area
for ii = par.Image0:par.dtImage:expTime
    fprintf('Time >>> %0.2f [min]\n',ii/60)
    
    
    for cc = 1:par.nw
        
        par.nfolder = dir([par.currentdirectory,'/Images/well_',num2str(cc),'/*.png']);%pictures taken with png format
        
        % reading current image
        if ii == par.Image0
            imNumb = 1;
        else
            imNumb = 1 + (ii - par.Image0)/15;
        end
        
        img_name = par.nfolder(imNumb).name;
        img_file = [par.currentdirectory,'/Images/well_',num2str(cc),'/',img_name];
        img_temp = imread(img_file);
        img_k{cc} = rgb2gray(img_temp);
        
        [dimY,dimX] = size(img_k{cc});
        
        % preparing mask
        areaMask = false(dimY,dimX);
        
        for yy = 1:dimY
            for xx = 1:dimX
                if img_k{cc}(yy,xx) > par.colorThreshold
                    areaMask(yy,xx) = true;
                end
            end
        end
        
        areaMask = imfill(areaMask,'holes');
        area(cc,imNumb) = sum(areaMask,'all');
    end
   
    
    %Display results
    g = figure(1);
    %set(gcf, 'Position', 0.9*get(0, 'Screensize'));
    axis tight manual %ensures that getframe() returns a consistent size
    clf;
    
    for well = 1:3
        tit = ['Well ',num2str(well)];
        
        % plotting image
        subplot(3,3,well), imshow(img_k{well}), title(tit)
        
        tgrid = 0:(par.dtMeasu/60):(par.dtMeasu/60)*((ii - 1));
        % plotting DP
        subplot(3,3,3 + well), 
                plot(tgrid,data{well}.deltaP(1:ii),'LineWidth',1.5,'MarkerSize',5);
                grid on
                
                %show one tick every 5 minutes 
                xticks(0:10:(par.dtMeasu/60)*(expTime - 1));
                xlim([0, (par.dtMeasu/60)*(expTime - 1)])
                ylim([0, 60])

                ylabel('DP [mbar]','FontSize',10)
                xlabel('time [min]','FontSize',10)

        % plotting Flowrate
        subplot(3,3,6 + well), 
                plot(tgrid,data{well}.flowrate(1:ii),'LineWidth',1.5,'MarkerSize',5);
                grid on
                
                %show one tick every 5 minutes 
                xticks(0:10:(par.dtMeasu/60)*(expTime - 1));
                xlim([0, (par.dtMeasu/60)*(expTime - 1)])
                ylim([4, 15])

                ylabel('Q [L/min]','FontSize',10)
                xlabel('time [min]','FontSize',10)

    end

    timeTag = ['time: ',num2str(ii/60,'%.2f'),' [min]'];
    dim = [.01 0.43 .3 .3];
    annotation('textbox',dim,'String',timeTag,'FitBoxToText','on');
	set(gcf,'color','w');

    %capturing frame
    frame = getframe(g);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    %write to the GIF file
    gifname = [filename,'.gif'];
    if ii == par.Image0
        imwrite(imind,cm,gifname,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,gifname,'gif','DelayTime',0.1,'WriteMode','append');
    end
    
    
end

savefig('FinalProfiles.fig')

%normalizing area
initArea = area(:,1);
areaN = area./initArea;

for cc = 1:par.nw
    data{cc}.areaPixel = area(cc,:);
    data{cc}.areaNPixel = areaN(cc,:);
end

%% Plotting

figure(2)

tgrid = 0:(par.dtImage/60):(par.dtImage/60)*(par.ni - 1);
label = {'Well 1','Well 2','Well 3'};

for well = 1:par.nw
    subplot(3,1,well)
    plot(tgrid,areaN(well,:),'k','Linewidth',1.5)
    
    grid on
    
    xticks(0:5:(par.ni - 1))
    yticks(0:0.25:1)
    
    ylim([0, 1.2])
    xlim([0, (par.dtImage/60)*par.ni])
    
    title(label{well})
    ylabel('Area [-]','FontSize',10)
    xlabel('time [min]','FontSize',10)
    
end

savefig('PixelCount.fig')

nameSave = filename;
save(nameSave,'data');

