clear variables; %close all;
set(0,'defaulttextinterpreter','latex'); rng('shuffle');
fontsize=6;
%% Parameters
figNum=0;
% general parameters
dt=0.0001;
% sf=dt; % sampling interval in s
% sf=0.0010;
% sf=0.0020;
% sf=0.0100;
sf=0.0200; % typical
HzSf=1; % window (and step) to calculate Hz in s
HzMax=160; % limit for input Hz axes
cbMax=1000; % scaling of color bar for imagesc plots
AMPAMax=1.5; % limit for AMPA axes
GABAMax=1.5; % limit for GABAR axes
availableNeurotransmitterMax=2.1;
directory='../../graphs/CorticoStriatal/';
% directory='../../../../../Data_iRIS/Data_actual/JamesH/crosstalk/0_5/';
fileExt='.dat';
% time frames
T=0:100;
Tperturbation=T; % same as above but for a perturbation
postprocess_InputSpikesFilter=true; % whether to filter input spikes
inputSpikeTimesT=90:T(end); % time window to filter
inputSpikeTimesTperturbation=0; % same as above but for a perturbation
postprocess_OutputSpikesFilter=true; % whether to filter output spikes
outputSpikeTimesT=inputSpikeTimesT; % time window to filter
outputSpikeTimesTperturbation=inputSpikeTimesTperturbation; % same as above but for a perturbation
finalMeasurementT=inputSpikeTimesT; % time window to calculate final measurements
finalMeasurementTperturbation=inputSpikeTimesTperturbation; % same as above but for a perturbation
% additional parameters
finalHeadroomStd=1; % accuracy of final headroom for determining time to headroom
headroomRisk=0.4; % the headroom demarcation line for risky synapses
%% Whether to load/process different data sources
postprocess_Glutamate = true;
postprocess_GABA = true;
glutamate_prep = 'Glutamate_';
GABA_prep = 'GABA_';

postprocess_PreIndexs = true;
postprocess_InputSpikes = true;
postprocess_Neurotransmitter = true;
postprocess_AvailableNeurotransmitter = true;
postprocess_CB1R = true;
postprocess_CB1Runbound = true;
postprocess_CB1Rcurrent = true;
postprocess_MAGL = true;
postprocess_GoodwinX = false;
postprocess_GoodwinY = false;
postprocess_GoodwinZ = false;
postprocess_GoodwinCannabinoidsX = true;
postprocess_GoodwinCannabinoidsk1 = true;
postprocess_CleftAstrocyteNeurotransmitter = true;
postprocess_CleftAstrocyteeCB = true;
postprocess_AMPAWeights = true;
postprocess_AMPA = true;
postprocess_mGluR5 = true;
postprocess_NMDARWeights = true;
postprocess_NMDAR = true;
postprocess_GABAWeights = true;
postprocess_GABAR = true;
postprocess_Ca = true;
postprocess_eCB = true;
postprocess_OutputSpikes = true;
postprocess_PostIndexs = true;

postprocess_Headroom = true;
postprocess_EIratio = false;%true;
postprocess_RiskySynapses = false;%true;

postprocess_Perturbation = false;%true;
postprocess_PerturbationHz = false;%true;
postprocess_PerturbationAMPA = false;%true;
postprocess_PerturbationGABA = false;%true;

Xdim = 500;%100;
Ydim = 1;
Zdim = 1;
glutamate_XdimInner = 25000; %40000;
glutamate_YdimInner = 1;
glutamate_ZdimInner = 1;
GABA_XdimInner = 25000;%10000;
GABA_YdimInner = 1;
GABA_ZdimInner = 1;
%%
for perturbation=0:(1*postprocess_Perturbation)
    clear glutamateInputSpike glutamateNeurotransmitter glutamateAvailableNeurotransmitter ...
        glutamatePreIndexs glutamatePreIndexs1D glutamateCB1R glutamateCB1Runbound ...
        glutamateCB1Rcurrent glutamateAMPA glutamateAMPAWeights ...
        glutamateAMPAWeights1D glutamatemGluR5 glutamateCa glutamateeCB ...
        GABAInputSpike GABANeurotransmitter GABAAvailableNeurotransmitter ...
        GABAPreIndexs GABAPreIndexs1D GABACB1R GABACB1Runbound ...
        GABACB1Rcurrent GABAMAGL GABAAMPA GABAAMPAWeights ...
        GABAAMPAWeights1D GABAmGluR5 GABACa GABAeCB ...
        outputSpike;
    % Change variables if a perturbation
    if (~perturbation)
        Trange=T;
        spikeTrange=inputSpikeTimesT;
        measurementTrange=finalMeasurementT;
    else
        Trange=Tperturbation;
        spikeTrange=inputSpikeTimesTperturbation;
        measurementTrange=finalMeasurementTperturbation;
    end
    sfRange = Trange(1)+sf:sf:Trange(end);
    % Load data
    if (postprocess_InputSpikes)
        if (postprocess_Glutamate)
            [glutamateInputSpike] = loadSpikes(directory, ...
                [glutamate_prep, 'PoissonSpikes'], fileExt, postprocess_InputSpikesFilter, ...
                spikeTrange(1), spikeTrange(end), Trange(1), Trange(end), dt, ...
                Xdim, Ydim, Zdim);
        end
        if (postprocess_GABA)
            [GABAInputSpike] = loadSpikes(directory, ...
                [GABA_prep, 'PoissonSpikes'], fileExt, postprocess_InputSpikesFilter, ...
                spikeTrange(1), spikeTrange(end), Trange(1), Trange(end), dt, ...
                Xdim, Ydim, Zdim);
        end
    end
    if (postprocess_Neurotransmitter) % do before PreIndexs to get XdimInner etc.
        if (postprocess_Glutamate)
            [glutamateNeurotransmitter] = ...
                load4D(directory, [glutamate_prep, 'BoutonNeurotransmitter'], ...
                fileExt, Trange(1), Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        end
        if (postprocess_GABA)
            [GABANeurotransmitter] = ...
                load4D(directory, [GABA_prep, 'BoutonNeurotransmitter'], ...
                fileExt, Trange(1), Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end 
    if (postprocess_PreIndexs) % return to correct order
        if (postprocess_Glutamate)
            [glutamatePreIndexs, glutamatePreIndexs1D] = load2D(directory, ...
                [glutamate_prep, 'BoutonIndexs'], fileExt, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        end
        if (postprocess_GABA)
            [GABAPreIndexs, GABAPreIndexs1D] = load2D(directory, ...
                [GABA_prep, 'BoutonIndexs'], fileExt, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end
    if (postprocess_AvailableNeurotransmitter)
        if (postprocess_Glutamate)
            [glutamateAvailableNeurotransmitter] = ...
                load4D(directory, [glutamate_prep, 'BoutonAvailableNeurotransmitter'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        end
        if (postprocess_GABA)
            [GABAAvailableNeurotransmitter] = ...
                load4D(directory, [GABA_prep, 'BoutonAvailableNeurotransmitter'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end 
    if (postprocess_CB1R)
        if (postprocess_Glutamate)
            [glutamateCB1R] = ...
                load4D(directory, [glutamate_prep, 'BoutonCB1R'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        end 
        if (postprocess_GABA)
            [GABACB1R] = ...
                load4D(directory, [GABA_prep, 'BoutonCB1R'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end 
    end
    if (postprocess_CB1Runbound)
        if (postprocess_Glutamate)
            [glutamateCB1Runbound] = ...
                load4D(directory, [glutamate_prep, 'BoutonCB1Runbound'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        end
        if (postprocess_GABA)
            [GABACB1Runbound] = ...
                load4D(directory, [GABA_prep, 'BoutonCB1Runbound'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end 
    if (postprocess_CB1Rcurrent)
        if (postprocess_Glutamate)
            [glutamateCB1Rcurrent] = ...
                load4D(directory, [glutamate_prep, 'BoutonCB1Rcurrent'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        end
        if (postprocess_GABA)
            [GABACB1Rcurrent] = ...
                load4D(directory, [GABA_prep, 'BoutonCB1Rcurrent'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end 
    if (postprocess_MAGL)
        if (postprocess_GABA)
            [GABAMAGL] = ...
                load4D(directory, [GABA_prep, 'BoutonMAGL'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end 
    if (postprocess_GoodwinX)
        if (postprocess_Glutamate)
            [glutamateGoodwinX] = ...
                load4D(directory, [glutamate_prep, 'Goodwin_X'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        end
        if (postprocess_GABA)
            [GABAGoodwinX] = ...
                load4D(directory, [GABA_prep, 'Goodwin_X'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end
    if (postprocess_GoodwinY)
        if (postprocess_Glutamate)
            [glutamateGoodwinY] = ...
                load4D(directory, [glutamate_prep, 'Goodwin_Y'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        end
        if (postprocess_GABA)
            [GABAGoodwinY] = ...
                load4D(directory, [GABA_prep, 'Goodwin_Y'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end
    if (postprocess_GoodwinZ)
        if (postprocess_Glutamate)
            [glutamateGoodwinZ] = ...
                load4D(directory, [glutamate_prep, 'Goodwin_Z'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        end
        if (postprocess_GABA)
            [GABAGoodwinZ] = ...
                load4D(directory, [GABA_prep, 'Goodwin_Z'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end   
    if (postprocess_GoodwinCannabinoidsX)
        if (postprocess_Glutamate)
            [glutamateGoodwinX] = ...
                load4D(directory, [glutamate_prep, 'GoodwinCannabinoids_X'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        end
        if (postprocess_GABA)
            [GABAGoodwinX] = ...
                load4D(directory, [GABA_prep, 'GoodwinCannabinoids_X'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end
    if (postprocess_GoodwinCannabinoidsk1)
        if (postprocess_GABA)
            [GABAGoodwink1] = ...
                load4D(directory, [GABA_prep, 'GoodwinCannabinoids_k1'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end 
    if (postprocess_CleftAstrocyteNeurotransmitter)
        if (postprocess_Glutamate)
            [glutamateCleftAstrocyteNeurotransmitter] = ...
                load4D(directory, ['CleftAstrocyteGlutamate'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        end
        if (postprocess_GABA)
            [GABACleftAstrocyteNeurotransmitter] = ...
                load4D(directory, ['CleftAstrocyteGABA'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end 
    if (postprocess_CleftAstrocyteeCB)
        if (postprocess_Glutamate)
            [glutamateCleftAstrocyteeCB] = ...
                load4D(directory, ['CleftAstrocyteGlutamateeCB'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        end
        if (postprocess_GABA)
            [GABACleftAstrocyteeCB] = ...
                load4D(directory, ['CleftAstrocyteGABAeCB'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end  
    if (postprocess_AMPAWeights && postprocess_Glutamate)
        if (perturbation && postprocess_PerturbationAMPA)
            [spineAMPAWeights, spineAMPAWeights1D] = load3D(directory, ...
                [glutamate_prep,'SpineAMPAweights_',num2str(Tperturbation(1)/dt)], fileExt, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        else
            [spineAMPAWeights, spineAMPAWeights1D] = load3D(directory, ...
                [glutamate_prep,'SpineAMPAweights_1'], fileExt, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        end
    end   
    if (postprocess_GABAWeights && postprocess_GABA)
        if (perturbation && postprocess_PerturbationGABA)
            [dendriteGABAWeights, dendriteGABAWeights1D] = load3D(directory, ...
                [GABA_prep,'DendriteGABAWeights_',num2str(Tperturbation(1)/dt)], fileExt, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        else
            [dendriteGABAWeights, dendriteGABAWeights1D] = load3D(directory, ...
                [GABA_prep,'DendriteGABAWeights_1'], fileExt, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end   
    if (postprocess_AMPA && postprocess_Glutamate)
        [spineAMPA] = ...
            load4D(directory, [glutamate_prep,'SpineAMPA'], fileExt, Trange(1), ...
            Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
    end 
    if (postprocess_GABAR && postprocess_GABA)
        [dendriteGABA] = ...
            load4D(directory, [GABA_prep,'DendriteGABA'], fileExt, Trange(1), ...
            Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
    end
    if (postprocess_mGluR5 && postprocess_Glutamate)
        if (~perturbation)
            [~, figNum] = newFigure(figNum, false);
            plotModulation(directory, [glutamate_prep,'SpinemGluR5modulation'], fileExt, ...
                'Spine mGluR5 Modulation Function', 'mGluR5', 'Ca2+', ...
                [glutamate_prep,'mGluR5modulation'], fontsize);
        end
        [spinemGluR5] = ...
            load4D(directory, [glutamate_prep,'SpinemGluR5'], fileExt, Trange(1), ...
            Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
    end
    if (postprocess_NMDARWeights && postprocess_Glutamate)
        [spineNMDARWeights, spineNMDARWeights1D] = load3D(directory, ...
            [glutamate_prep,'SpineNMDARweights_1'], fileExt, ...
            glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
    end   
    if (postprocess_NMDAR && postprocess_Glutamate)
        [spineNMDAR] = ...
            load4D(directory, [glutamate_prep,'SpineNMDARcurrent'], fileExt, Trange(1), ...
            Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
    end
    if (postprocess_Ca)
        if (postprocess_Glutamate)
            [spineCa] = ...
                load4D(directory, [glutamate_prep,'SpineCa'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        end
        if (postprocess_GABA)
            [dendriteCa] = ...
                load4D(directory, [GABA_prep,'DendriteCa'], fileExt, Trange(1), ...
                Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
        end
    end
    if (postprocess_eCB)
        if (~perturbation)
            if (postprocess_Glutamate)
                [~, figNum] = newFigure(figNum, false);
                plotModulation(directory, [glutamate_prep,'SpineeCBproduction'], fileExt, ...
                    'Spine eCB Modulation Function', 'Ca2+', 'eCB+', ...
                    [glutamate_prep,'eCBproduction'], fontsize);
            end
            if (postprocess_GABA)
                [~, figNum] = newFigure(figNum, false);
                plotModulation(directory, [GABA_prep,'DendriteeCBproduction'], fileExt, ...
                    'Dendrite eCB Modulation Function', 'Ca2+', 'eCB+', ...
                    [GABA_prep,'eCBproduction'], fontsize);
            end
        end
        [spineeCB] = ...
            load4D(directory, [glutamate_prep,'SpineeCB'], fileExt, Trange(1), ...
            Trange(end), sf, ...
                glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner);
        [dendriteeCB] = ...
            load4D(directory, [GABA_prep,'DendriteeCB'], fileExt, Trange(1), ...
            Trange(end), sf, ...
                GABA_XdimInner, GABA_YdimInner, GABA_ZdimInner);
    end
    if (postprocess_OutputSpikes)
        [outputSpike] = loadSpikes(directory, ...
            'Output_Spikes', fileExt, postprocess_OutputSpikesFilter, ...
            spikeTrange(1), spikeTrange(end), Trange(1), Trange(end), dt, ...
                Xdim, Ydim, Zdim);
    end
    if (postprocess_PostIndexs)
        if (postprocess_Glutamate)
            fid = fopen([directory,'OutputAMPAIndexs',fileExt],'r');
            i = 1;
            OutputAMPAIndexs1D = [];
            for x=1:Xdim
                for y=1:Ydim
                    for z=1:Zdim
                        OutputAMPAIndexs_N(i) = fread(fid, 1, 'int');
                        OutputAMPAIndexs{i} = fread(fid, OutputAMPAIndexs_N(i)*2, 'int');
                        OutputAMPAIndexs{i}(1:2:end) = [];
                        OutputAMPAIndexs1D = [OutputAMPAIndexs1D; OutputAMPAIndexs{i}(:)];
                        i = i + 1;
                    end
                end
            end
            fclose(fid);
            clear fid; 
        end
        if (postprocess_GABA)
            fid = fopen([directory,'OutputGABARIndexs',fileExt],'r');
            i = 1;
            OutputGABARIndexs1D = [];
            for x=1:Xdim
                for y=1:Ydim
                    for z=1:Zdim
                        OutputGABARIndexs_N(i) = fread(fid, 1, 'int');
                        OutputGABARIndexs{i} = fread(fid, OutputGABARIndexs_N(i)*2, 'int');
                        OutputGABARIndexs{i}(1:2:end) = [];
                        OutputGABARIndexs1D = [OutputGABARIndexs1D; OutputGABARIndexs{i}(:)];
                        i = i + 1;
                    end
                end
            end
            fclose(fid);
            clear fid; 
        end
    end
    % Plot
    if (postprocess_InputSpikes)
        if (postprocess_Glutamate)
            [~, figNum] = newFigure(figNum, false); % Population firing rate histogram over all time
            glutamateInHz = calculateHz(glutamateInputSpike, spikeTrange(1), ...
                spikeTrange(end), HzSf, Xdim, ...
                Ydim, Zdim, dt);
            histogram(mean(glutamateInHz(:),2));
            title(['Glutamate in Hz',perturbationString(perturbation,0,1)]); 
            xlabel('Hz'); ylabel('count');
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            print([directory,glutamate_prep,'inHz',perturbationString(perturbation,1,0)],'-dpng');
            print([directory,glutamate_prep,'inHz',perturbationString(perturbation,1,0)],'-depsc');
        end
        if (postprocess_GABA)
            [~, figNum] = newFigure(figNum, false);
            GABAInHz = calculateHz(GABAInputSpike, spikeTrange(1), ...
                spikeTrange(end), HzSf, Xdim, ...
                Ydim, Zdim, dt);
            histogram(mean(GABAInHz(:),2));
            title(['GABA in Hz',perturbationString(perturbation,0,1)]); 
            xlabel('Hz'); ylabel('count');
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            print([directory,GABA_prep,'inHz',perturbationString(perturbation,1,0)],'-dpng');
            print([directory,GABA_prep,'inHz',perturbationString(perturbation,1,0)],'-depsc');
        end
    end
    if (glutamate_XdimInner < GABA_XdimInner)
        synapse=randi(glutamate_XdimInner,1,1); % only for Xdimension
    else
        synapse=randi(GABA_XdimInner,1,1); % only for Xdimension
    end
    [~, figNum] = newFigure(figNum, false);
    glutamateFigNum = figNum;
    [~, figNum] = newFigure(figNum, false);
    GABAFigNum = figNum;
    if (postprocess_InputSpikes)
        if (postprocess_Glutamate)
            figure(glutamateFigNum);
            subplot(4,5,1); 
            scatter(glutamateInputSpike{glutamatePreIndexs(1,synapse,1,1),1,1}.*dt, ...
                ones(numel(glutamateInputSpike{glutamatePreIndexs(1,synapse,1,1),1,1}),1),'.');
            title(['Glutamate Input Spikes',perturbationString(perturbation,0,1)]); 
            xlim([Trange(1) Trange(end)]);
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,1); 
            scatter(GABAInputSpike{GABAPreIndexs(1,synapse,1,1),1,1}.*dt, ...
                ones(numel(GABAInputSpike{GABAPreIndexs(1,synapse,1,1),1,1}),1),'.');
            title(['GABA Input Spikes',perturbationString(perturbation,0,1)]); 
            xlim([Trange(1) Trange(end)]);
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
    end
    if (postprocess_Neurotransmitter || postprocess_AvailableNeurotransmitter)
        if (postprocess_Glutamate)
            figure(glutamateFigNum);
            subplot(4,5,2); hold on;
            if (postprocess_Neurotransmitter)
                plot(sfRange, glutamateNeurotransmitter(:,synapse,1,1));
            end
            if (postprocess_AvailableNeurotransmitter)
                plot(sfRange, glutamateAvailableNeurotransmitter(:,synapse,1,1));
            end
            xlim([Trange(1) Trange(end)]);
            title(['Glutamate and Available Glutamate',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,2); hold on;
            if (postprocess_Neurotransmitter)
                plot(sfRange, GABANeurotransmitter(:,synapse,1,1));
            end
            if (postprocess_AvailableNeurotransmitter)
                plot(sfRange, GABAAvailableNeurotransmitter(:,synapse,1,1));
            end
            xlim([Trange(1) Trange(end)]);
            title(['GABA and Available GABA',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
    end
    if (postprocess_CB1R)
        if (postprocess_Glutamate)
            figure(glutamateFigNum);
            subplot(4,5,3);
            plot(sfRange, glutamateCB1R(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['Glutamate CB1R',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,3);
            plot(sfRange, GABACB1R(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['GABA CB1R',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
    end 
    if (postprocess_CB1Runbound)
        if (postprocess_Glutamate)
            figure(glutamateFigNum);
            subplot(4,5,4);
            plot(sfRange, glutamateCB1Runbound(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['Glutamate CB1R unbound',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,4);
            plot(sfRange, GABACB1Runbound(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['GABA CB1R unbound',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
    end 
    if (postprocess_CB1Rcurrent)
        if (postprocess_Glutamate)
            figure(glutamateFigNum);
            subplot(4,5,5);
            plot(sfRange, glutamateCB1Rcurrent(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['Glutamate CB1R Current',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,5);
            plot(sfRange, GABACB1Rcurrent(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['GABA CB1R Current',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
    end 
    if (postprocess_MAGL)
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,6);
            plot(sfRange, GABAMAGL(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['GABA MAGL',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
    end 
    if (postprocess_GoodwinX || postprocess_GoodwinCannabinoidsX)
        if (postprocess_Glutamate)
            figure(glutamateFigNum);
            subplot(4,5,7);
            plot(sfRange, glutamateGoodwinX(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['Glutamate Goodwin X',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,7);
            plot(sfRange, GABAGoodwinX(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['GABA Goodwin X',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
    end 
    if (postprocess_GoodwinY)
        if (postprocess_Glutamate)
            figure(glutamateFigNum);
            subplot(4,5,8);
            plot(sfRange, glutamateGoodwinY(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['Glutamate Goodwin Y',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,8);
            plot(sfRange, GABAGoodwinY(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['GABA Goodwin Y',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
    end 
    if (postprocess_GoodwinZ)
        if (postprocess_Glutamate)
            figure(glutamateFigNum);
            subplot(4,5,9);
            plot(sfRange, glutamateGoodwinZ(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['Glutamate Goodwin Z',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end 
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,9);
            plot(sfRange, GABAGoodwinZ(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['GABA Goodwin Z',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end 
    end
    if (postprocess_GoodwinCannabinoidsk1)
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,10);
            plot(sfRange, GABAGoodwink1(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['GABA Goodwin k1',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end 
    end
    if (postprocess_CleftAstrocyteNeurotransmitter)
        if (postprocess_Glutamate)
            figure(glutamateFigNum);
            subplot(4,5,11);
            plot(sfRange, glutamateCleftAstrocyteNeurotransmitter(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['Glutamate Cleft',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,11);
            plot(sfRange, GABACleftAstrocyteNeurotransmitter(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['GABA Cleft',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
    end
    if (postprocess_CleftAstrocyteeCB)
        if (postprocess_Glutamate)
            figure(glutamateFigNum);
            subplot(4,5,12);
            plot(sfRange, glutamateCleftAstrocyteeCB(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['Glutamate Cleft eCB',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,12);
            plot(sfRange, GABACleftAstrocyteeCB(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['GABA Cleft eCB',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
    end    
    if ((postprocess_AMPAWeights || postprocess_AMPA) && postprocess_Glutamate)
        figure(glutamateFigNum);
        subplot(4,5,13); hold on;
        if (postprocess_AMPAWeights)
            plot(sfRange, ones(1,numel(sfRange))*spineAMPAWeights(synapse,1,1));
            yyaxis right;
        end
        if (postprocess_AMPA)
            plot(sfRange, spineAMPA(:,synapse,1,1));  
        end
        xlim([Trange(1) Trange(end)]);
        title(['Spine AMPA/AMPAWeights',perturbationString(perturbation,0,1)]); 
        set(gca,'XTick',[]);
        set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    end 
    if ((postprocess_GABAWeights || postprocess_GABAR) && postprocess_GABA)
        figure(GABAFigNum);
        subplot(4,5,13); hold on;
        if (postprocess_GABAWeights)
            plot(sfRange, ones(1,numel(sfRange))*dendriteGABAWeights(synapse,1,1));
            yyaxis right;
        end
        if (postprocess_GABAR)
            plot(sfRange, dendriteGABA(:,synapse,1,1));  
        end
        xlim([Trange(1) Trange(end)]);
        title(['Dendrite GABA/GABAWeights',perturbationString(perturbation,0,1)]); 
        set(gca,'XTick',[]);
        set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    end
    if (postprocess_mGluR5 && postprocess_Glutamate)
        figure(glutamateFigNum);
        subplot(4,5,14);
        plot(sfRange, spinemGluR5(:,synapse,1,1));
        xlim([Trange(1) Trange(end)]);
        title(['Spine mGluR5',perturbationString(perturbation,0,1)]); 
        set(gca,'XTick',[]);
        set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    end
    if ((postprocess_NMDARWeights || postprocess_NMDAR) && postprocess_Glutamate)
        figure(glutamateFigNum);
        subplot(4,5,15); hold on;
        if (postprocess_NMDARWeights)
            plot(sfRange, ones(1,numel(sfRange))*spineNMDARWeights(synapse,1,1));
            yyaxis right;            
        end
        if (postprocess_NMDAR)
            plot(sfRange, spineNMDAR(:,synapse,1,1));
        end
        xlim([Trange(1) Trange(end)]);
        title(['Spine NMDAR/NMDARWeights',perturbationString(perturbation,0,1)]);
        set(gca,'XTick',[]);
        set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    end
    if (postprocess_Ca)
        if (postprocess_Glutamate)
            figure(glutamateFigNum);
            subplot(4,5,16);
            plot(sfRange, spineCa(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['Spine Ca',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,16);
            plot(sfRange, dendriteCa(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['Dendrite Ca',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
    end
    if (postprocess_eCB)
        if (postprocess_Glutamate)
            figure(glutamateFigNum);
            subplot(4,5,17);
            plot(sfRange, spineeCB(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['Spine eCB',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,17);
            plot(sfRange, dendriteeCB(:,synapse,1,1));
            xlim([Trange(1) Trange(end)]);
            title(['Dendrite eCB',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
    end
    if (postprocess_OutputSpikes)
        if (postprocess_Glutamate)
            figure(glutamateFigNum);
            subplot(4,5,18);
            scatter(outputSpike{glutamatePreIndexs(1,synapse,1,1),1,1}.*dt, ...
                ones(numel(outputSpike{glutamatePreIndexs(1,synapse,1,1),1,1}),1),'.');
            xlim([Trange(1) Trange(end)]);
            title(['Glutamate Output Spike',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
        if (postprocess_GABA)
            figure(GABAFigNum);
            subplot(4,5,18);
            scatter(outputSpike{GABAPreIndexs(1,synapse,1,1),1,1}.*dt, ...
                ones(numel(outputSpike{GABAPreIndexs(1,synapse,1,1),1,1}),1),'.');
            xlim([Trange(1) Trange(end)]);
            title(['GABA Output Spike',perturbationString(perturbation,0,1)]); 
            set(gca,'XTick',[]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        end
    end
    if (postprocess_InputSpikes || postprocess_Neurotransmitter ...
            || postprocess_AvailableNeurotransmitter || postprocess_CB1R ...
            || postprocess_CB1Runbound || postprocess_CB1Rcurrent ...
            || postprocess_CleftNeurotransmitter || postprocess_AMPA ...
            || postprocess_AMPAWeights || postprocess_GABAR ...
            || postprocess_GABAWeights || postprocess_mGluR5 ...
            || postprocess_Ca || postprocess_eCB ...
            || postprocess_OutputSpikes)
        if (postprocess_Glutamate)
            figure(glutamateFigNum);
            print([directory,glutamate_prep,'Components',perturbationString(perturbation,1,0)],'-dpng');
            print([directory,glutamate_prep,'Components',perturbationString(perturbation,1,0)],'-depsc');
        end
        if (postprocess_GABA)
            figure(GABAFigNum);
            print([directory,GABA_prep,'Components',perturbationString(perturbation,1,0)],'-dpng');
            print([directory,GABA_prep,'Components',perturbationString(perturbation,1,0)],'-depsc');
        end
    end
    if (postprocess_OutputSpikes)
        [~, figNum] = newFigure(figNum, false); % Population firing rate histogram over all time
        outHz = calculateHz(outputSpike, spikeTrange(1), ...
            spikeTrange(end), HzSf, Xdim, ...
            Ydim, Zdim, dt);
        histogram(mean(outHz(:),2));
        title(['Out population Hz',perturbationString(perturbation,0,1)]); 
        xlabel('Hz'); ylabel('count');
        set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        print([directory,'outPopulationHz',perturbationString(perturbation,1,0)],'-dpng');
        print([directory,'outPopulationHz',perturbationString(perturbation,1,0)],'-depsc');
        
        [~, figNum] = newFigure(figNum, false); % Individual firing rates
        plot(outHz);
        title(['Out individual Hz',perturbationString(perturbation,0,1)]); 
        xlabel('time'); ylabel('Hz');
        set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        print([directory,'outIndividualHz',perturbationString(perturbation,1,0)],'-dpng');
        print([directory,'outIndividualHz',perturbationString(perturbation,1,0)],'-depsc');
        
%         [~, figNum] = newFigure(figNum, false); % Mean individual firing rate
%         outHzMean = mean(outHz,2);
%         outHzStd = std(outHz,0,2);
%         X=[outputSpikeTimesT(1):HzSf:outputSpikeTimesT(end),fliplr(outputSpikeTimesT(1):HzSf:outputSpikeTimesT(end))];
%         Y=[outHzMean'+outHzStd',fliplr(outHzMean'-outHzStd')];
%         fill(X,Y,'b','LineStyle','none');
%         hold on;
%         plot(outputSpikeTimesT(1):HzSf:outputSpikeTimesT(end),outHzMean,'b','LineWidth',1.5);
%         alpha(0.35);
%         title(['Out individual mean Hz',perturbationString(perturbation,0,1)]); 
%         xlabel('time'); ylabel('Hz');
%         set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
%         print([directory,'outIndividualMeanHz',perturbationString(perturbation,1,0)],'-dpng');
%         print([directory,'outIndividualMeanHz',perturbationString(perturbation,1,0)],'-depsc');
    end    
    if (postprocess_AMPA && postprocess_NMDAR)
        [~, figNum] = newFigure(figNum, false);
        
        subplot(1,2,1);
        histogram(spineAMPA);
        title('AMPA');
        xlim([0 0.02]);
        subplot(1,2,2);
        histogram(spineNMDAR);
        title('NMDA');
        xlim([0 0.02]);
        
        print([directory,glutamate_prep,'AMPA_NMDAR',perturbationString(perturbation,1,0)],'-dpng');
        print([directory,glutamate_prep,'AMPA_NMDAR',perturbationString(perturbation,1,0)],'-depsc');
    end 
    if (postprocess_mGluR5 && postprocess_Ca && postprocess_eCB && postprocess_Glutamate)
        [~, figNum] = newFigure(figNum, false);
        
        Hrange = 0:max(spinemGluR5(:))/100:max(spinemGluR5(:));
        spinemGluR5Hist = hist3D(Hrange, spinemGluR5);
        subplot(3,1,1);
        imagesc([0 size(spinemGluR5,1)],Hrange,spinemGluR5Hist);
        title(['Spine mGluR5',perturbationString(perturbation,0,1)]); 
        xlabel('time'); ylabel('mGluR5');
        colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
        set(gca,'ydir','normal');
        set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);        

        Hrange = 0:max(spineCa(:))/100:max(spineCa(:));
        spineCaHist = hist3D(Hrange, spineCa);
        subplot(3,1,2);
        imagesc([0 size(spineCa,1)],Hrange,spineCaHist);
        title(['Spine Ca2+',perturbationString(perturbation,0,1)]);
        xlabel('time'); ylabel('Ca2+');
        colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
        set(gca,'ydir','normal');
        set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);

        Hrange = 0:0.01:1;
        spineeCBHist = hist3D(Hrange, spineeCB);
        subplot(3,1,3);
        imagesc([0 size(spineeCB,1)],Hrange,spineeCBHist);
        title(['Spine eCB',perturbationString(perturbation,0,1)]);
        xlabel('time'); ylabel('eCB'); ylim([0 1]);
        colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
        set(gca,'ydir','normal');
        set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        
        print([directory,glutamate_prep,'spine_mGluR5_Ca_eCB',perturbationString(perturbation,1,0)],'-dpng');
        print([directory,glutamate_prep,'spine_mGluR5_Ca_eCB',perturbationString(perturbation,1,0)],'-depsc');
        clear Hrange spinemGluR5Hist spineCaHist spineeCBHist;
    end    
    if (postprocess_Ca && postprocess_eCB && postprocess_GABA)
        [~, figNum] = newFigure(figNum, false);

        Hrange = 0:max(dendriteCa(:))/100:max(dendriteCa(:));
        dendriteCaHist = hist3D(Hrange, dendriteCa);
        subplot(2,1,1);
        imagesc([0 size(dendriteCa,1)],Hrange,dendriteCaHist);
        title(['Dendrite Ca2+',perturbationString(perturbation,0,1)]);
        xlabel('time'); ylabel('Ca2+');
        colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
        set(gca,'ydir','normal');
        set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);

        Hrange = 0:0.01:1;
        dendriteeCBHist = hist3D(Hrange, dendriteeCB);
        subplot(2,1,2);
        imagesc([0 size(dendriteeCB,1)],Hrange,dendriteeCBHist);
        title(['Dendrite eCB',perturbationString(perturbation,0,1)]);
        xlabel('time'); ylabel('eCB'); ylim([0 1]);
        colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
        set(gca,'ydir','normal');
        set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
        
        print([directory,GABA_prep,'dendrite_Ca_eCB',perturbationString(perturbation,1,0)],'-dpng');
        print([directory,GABA_prep,'dendrite_Ca_eCB',perturbationString(perturbation,1,0)],'-depsc');
        clear Hrange dendriteCaHist dendriteeCBHist;
    end 
    if (postprocess_CB1R || postprocess_CB1Runbound || postprocess_MAGL)
        if (postprocess_Glutamate)
            [~, figNum] = newFigure(figNum, false);
            Hrange = 0:1/50:1;

            if (postprocess_CB1R)
                glutamateCB1RHist = hist3D(Hrange, glutamateCB1R);
                subplot(2,1,1);
                imagesc([0 size(glutamateCB1R,1)],Hrange,glutamateCB1RHist);
                title(['Glutamate CB1R',perturbationString(perturbation,0,1)]); 
                xlabel('time'); ylabel('CB1R');
                colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
                set(gca,'ydir','normal');
                set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            end

            if (postprocess_CB1Runbound)
                glutamateCB1RunboundHist = hist3D(Hrange, glutamateCB1Runbound);
                subplot(2,1,2);
                imagesc([0 size(glutamateCB1Runbound,1)],Hrange,glutamateCB1RunboundHist);
                title(['Glutamate Unbound CB1R',perturbationString(perturbation,0,1)]); 
                xlabel('time'); ylabel('Unbound CB1R');
                colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
                set(gca,'ydir','normal');
                set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            end
            
            print([directory,glutamate_prep,'CB1R_CB1Runbound',perturbationString(perturbation,1,0)],'-dpng');
            print([directory,glutamate_prep,'CB1R_CB1Runbound',perturbationString(perturbation,1,0)],'-depsc');
            clear Hrange glutamateCB1RHist glutamateCB1RunboundHist; 
        end
        if (postprocess_GABA)
            [~, figNum] = newFigure(figNum, false);
            Hrange = 0:1/50:1;

            if (postprocess_CB1R)
                GABACB1RHist = hist3D(Hrange, GABACB1R);
                subplot(3,1,1);
                imagesc([0 size(GABACB1R,1)],Hrange,GABACB1RHist);
                title(['GABA CB1R',perturbationString(perturbation,0,1)]); 
                xlabel('time'); ylabel('CB1R');
                colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
                set(gca,'ydir','normal');
                set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            end

            if (postprocess_CB1Runbound)
                GABACB1RunboundHist = hist3D(Hrange, GABACB1Runbound);
                subplot(3,1,2);
                imagesc([0 size(GABACB1Runbound,1)],Hrange,GABACB1RunboundHist);
                title(['GABA Unbound CB1R',perturbationString(perturbation,0,1)]); 
                xlabel('time'); ylabel('Unbound CB1R');
                colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
                set(gca,'ydir','normal');
                set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            end

            if (postprocess_MAGL)
                MAGLHist = hist3D(Hrange, GABAMAGL);
                subplot(3,1,3);
                imagesc([0 size(GABAMAGL,1)],Hrange,MAGLHist);
                title(['GABA MAGL',perturbationString(perturbation,0,1)]); 
                xlabel('time'); ylabel('MAGL');
                colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
                set(gca,'ydir','normal');
                set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            end
            
            print([directory,GABA_prep,'CB1R_CB1Runbound_MAGL',perturbationString(perturbation,1,0)],'-dpng');
            print([directory,GABA_prep,'CB1R_CB1Runbound_MAGL',perturbationString(perturbation,1,0)],'-depsc');
            clear Hrange GABACB1RHist GABACB1RunboundHist MAGLHist; 
        end
    end
    if (postprocess_GoodwinX || postprocess_GoodwinCannabinoidsX || postprocess_GoodwinCannabinoidsk1)
        if (postprocess_Glutamate)
            [~, figNum] = newFigure(figNum, false);
            Hrange = 0:10/50:10;

            if (postprocess_GoodwinX || postprocess_GoodwinCannabinoidsX)
                glutamateGoodwinXHist = hist3D(Hrange, glutamateGoodwinX);
                imagesc([0 size(glutamateGoodwinX,1)],Hrange,glutamateGoodwinXHist);
                title(['Glutamate Goodwin X',perturbationString(perturbation,0,1)]); 
                xlabel('time'); ylabel('X');
                colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
                set(gca,'ydir','normal');
                set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            end
            
            print([directory,glutamate_prep,'GoodwinX_Goodwink1',perturbationString(perturbation,1,0)],'-dpng');
            print([directory,glutamate_prep,'GoodwinX_Goodwink1',perturbationString(perturbation,1,0)],'-depsc');
            clear Hrange glutamateGoodwinXHist; 
        end
        if (postprocess_GABA)
            [~, figNum] = newFigure(figNum, false);
            Hrange = 0:10/50:10;

            if (postprocess_GoodwinX || postprocess_GoodwinCannabinoidsX)
                GABAGoodwinXHist = hist3D(Hrange, GABAGoodwinX);
                subplot(2,1,1);
                imagesc([0 size(GABAGoodwinX,1)],Hrange,GABAGoodwinXHist);
                title(['GABA Goodwin X',perturbationString(perturbation,0,1)]); 
                xlabel('time'); ylabel('X');
                colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
                set(gca,'ydir','normal');
                set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            end

            if (postprocess_GoodwinCannabinoidsk1)
                GABAGoodwink1Hist = hist3D(Hrange, GABAGoodwink1);
                subplot(2,1,2);
                imagesc([0 size(GABAGoodwink1,1)],Hrange,GABAGoodwink1Hist);
                title(['GABA Goodwin k1',perturbationString(perturbation,0,1)]); 
                xlabel('time'); ylabel('k1');
                colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
                set(gca,'ydir','normal');
                set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            end
            
            print([directory,GABA_prep,'GoodwinX_Goodwink1',perturbationString(perturbation,1,0)],'-dpng');
            print([directory,GABA_prep,'GoodwinX_Goodwink1',perturbationString(perturbation,1,0)],'-depsc');
            clear Hrange GABAGoodwinXHist GABAGoodwink1Hist; 
        end
    end
    if (postprocess_AvailableNeurotransmitter)
        if (postprocess_Glutamate)
            Hrange = 0:availableNeurotransmitterMax/100:availableNeurotransmitterMax;
            glutamateAvailableNeurotransmitterHist = hist3D(Hrange, ...
                glutamateAvailableNeurotransmitter);
            [~, figNum] = newFigure(figNum, false);
            imagesc([0 size(glutamateAvailableNeurotransmitter,1)]...
                ,Hrange,glutamateAvailableNeurotransmitterHist);
            title(['Glutamate Available Neurotransmitter',perturbationString(perturbation,0,1)]);
            xlabel('time'); 
            ylabel('Available Neurotransmitter');
            colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
            set(gca,'ydir','normal');
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            print([directory,glutamate_prep,'availableNeurotransmitter_hist',perturbationString(perturbation,1,0)],'-dpng');
            print([directory,glutamate_prep,'availableNeurotransmitter_hist',perturbationString(perturbation,1,0)],'-depsc');
            clear Hrange glutamateAvailableNeurotransmitterHist;
        end
        if (postprocess_GABA)
            Hrange = 0:availableNeurotransmitterMax/100:availableNeurotransmitterMax;
            GABAAvailableNeurotransmitterHist = hist3D(Hrange, ...
                GABAAvailableNeurotransmitter);
            [~, figNum] = newFigure(figNum, false);
            imagesc([0 size(GABAAvailableNeurotransmitter,1)]...
                ,Hrange,GABAAvailableNeurotransmitterHist);
            title(['GABA Available Neurotransmitter',perturbationString(perturbation,0,1)]);
            xlabel('time'); 
            ylabel('Available Neurotransmitter');
            colormap(flipud(hot)); colorbar(); caxis([0 cbMax]);
            set(gca,'ydir','normal');
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            print([directory,GABA_prep,'availableNeurotransmitter_hist',perturbationString(perturbation,1,0)],'-dpng');
            print([directory,GABA_prep,'availableNeurotransmitter_hist',perturbationString(perturbation,1,0)],'-depsc');
            clear Hrange GABAAvailableNeurotransmitterHist;
        end
    end
    if (postprocess_Headroom)
        if (postprocess_Glutamate)
            Hrange = -availableNeurotransmitterMax:availableNeurotransmitterMax/200:availableNeurotransmitterMax;
            [~, figNum] = newFigure(figNum, false);
            glutamateHeadroomH = histHeadroom3D(Hrange, glutamateAvailableNeurotransmitter, ...
                spineAMPAWeights1D);
            imagesc([0 size(glutamateAvailableNeurotransmitter,1)],Hrange,glutamateHeadroomH);
            title(['Glutamate Excess Neurotransmitter',perturbationString(perturbation,0,1)]);
            xlabel('time'); ylabel('excess neurotransmitter');
            colormap(flipud(hot)); colorbar();
            set(gca,'ydir','normal');
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            print([directory,glutamate_prep,'excessNeurotransmitter',...
                perturbationString(perturbation,1,0)],'-dpng');
            print([directory,glutamate_prep,'excessNeurotransmitter',...
                perturbationString(perturbation,1,0)],'-depsc');
            [~, figNum] = newFigure(figNum, false);
            plot(Hrange,glutamateHeadroomH(:,1));
            title(['Glutamate Excess Neurotransmitter Distribution - Initial Condition',...
                perturbationString(perturbation,0,1)]);
            xlabel('excess neurotransmitter'); ylabel('count');
            xlim([Hrange(1) Hrange(end)]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            print([directory,glutamate_prep,'excessNeurotransmitterDistInitialCondition'...
                ,perturbationString(perturbation,1,0)],'-dpng');
            print([directory,glutamate_prep,'excessNeurotransmitterDistInitialCondition'...
                ,perturbationString(perturbation,1,0)],'-depsc');
            [~, figNum] = newFigure(figNum, false);
            plot(Hrange,glutamateHeadroomH(:,end));
            title(['Glutamate Excess Neurotransmitter Distribution - Steady State',...
                perturbationString(perturbation,0,1)]);
            xlabel('excess neurotransmitter'); ylabel('count');
            xlim([Hrange(1) Hrange(end)]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            print([directory,glutamate_prep,'excessNeurotransmitterDistSteadyState',...
                perturbationString(perturbation,1,0)],'-dpng');
            print([directory,glutamate_prep,'excessNeurotransmitterDistSteadyState',...
                perturbationString(perturbation,1,0)],'-depsc');
            clear Hrange glutamateHeadroomH;
        end
        if (postprocess_GABA)
            Hrange = -availableNeurotransmitterMax:availableNeurotransmitterMax/200:availableNeurotransmitterMax;
            [~, figNum] = newFigure(figNum, false);
            GABAHeadroomH = histHeadroom3D(Hrange, GABAAvailableNeurotransmitter, ...
                dendriteGABAWeights1D);
            imagesc([0 size(GABAAvailableNeurotransmitter,1)],Hrange,GABAHeadroomH);
            title(['GABA Excess Neurotransmitter',perturbationString(perturbation,0,1)]);
            xlabel('time'); ylabel('excess neurotransmitter');
            colormap(flipud(hot)); colorbar();
            set(gca,'ydir','normal');
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            print([directory,GABA_prep,'excessNeurotransmitter',...
                perturbationString(perturbation,1,0)],'-dpng');
            print([directory,GABA_prep,'excessNeurotransmitter',...
                perturbationString(perturbation,1,0)],'-depsc');
            [~, figNum] = newFigure(figNum, false);
            plot(Hrange,GABAHeadroomH(:,1));
            title(['GABA Excess Neurotransmitter Distribution - Initial Condition',...
                perturbationString(perturbation,0,1)]);
            xlabel('excess neurotransmitter'); ylabel('count');
            xlim([Hrange(1) Hrange(end)]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            print([directory,GABA_prep,'excessNeurotransmitterDistInitialCondition'...
                ,perturbationString(perturbation,1,0)],'-dpng');
            print([directory,GABA_prep,'excessNeurotransmitterDistInitialCondition'...
                ,perturbationString(perturbation,1,0)],'-depsc');
            [~, figNum] = newFigure(figNum, false);
            plot(Hrange,GABAHeadroomH(:,end));
            title(['GABA Excess Neurotransmitter Distribution - Steady State',...
                perturbationString(perturbation,0,1)]);
            xlabel('excess neurotransmitter'); ylabel('count');
            xlim([Hrange(1) Hrange(end)]);
            set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
            print([directory,GABA_prep,'excessNeurotransmitterDistSteadyState',...
                perturbationString(perturbation,1,0)],'-dpng');
            print([directory,GABA_prep,'excessNeurotransmitterDistSteadyState',...
                perturbationString(perturbation,1,0)],'-depsc');
            clear Hrange GABAHeadroomH;
        end
        if (postprocess_Glutamate && postprocess_GABA)
            [...
            glutamateHeadroomFinal, glutamateHeadroomFinalStd, ...
            glutamateHzFinal, glutamateTimeToHeadroom, glutamateTimeToHeadroomConverged, ...
            glutamateCB1RFinal, glutamateAvailableNeurotransmitterFinal, spineCaFinal, ...
            glutamateCleftAstrocyteeCBFinal, ...
            ...
            GABAHeadroomFinal, GABAHeadroomFinalStd, ...
            GABAHzFinal, GABATimeToHeadroom, GABATimeToHeadroomConverged, ...
            GABACB1RFinal, GABAAvailableNeurotransmitterFinal, dendriteCaFinal, ...
            GABACleftAstrocyteeCBFinal, ...
            ...
            outHzFinal, ...
            ...
            figNum] = ...
            plotComponentVsHeadroom(figNum, Xdim, Ydim, Zdim, ...
            glutamate_XdimInner, glutamate_YdimInner, glutamate_ZdimInner, ...
            ...
            glutamateAvailableNeurotransmitter, spineAMPAWeights, spineAMPAWeights1D, glutamateInHz, ...
            glutamatePreIndexs1D, AMPAMax, glutamateCB1R, spineCa, ...
            glutamateCleftAstrocyteeCB, OutputAMPAIndexs1D, ...
            ...
            GABAAvailableNeurotransmitter, dendriteGABAWeights, dendriteGABAWeights1D, GABAInHz, ...
            GABAPreIndexs1D, GABAMax, GABACB1R, dendriteCa, ...
            GABACleftAstrocyteeCB, OutputGABARIndexs1D, ...
            ...
            outHz, ...
            ...
            measurementTrange, measurementTrange, finalHeadroomStd, Trange(1), sf, HzSf, ...
            spikeTrange(1), HzMax, ...
            directory, glutamate_prep, GABA_prep, fontsize);
        end
        if (postprocess_RiskySynapses)
            if (postprocess_Glutamate)
                [~, figNum] = newFigure(figNum, true);
                h = scatterhist(glutamateHeadroomFinal(glutamateTimeToHeadroomConverged), ...
                    glutamateTimeToHeadroom(glutamateTimeToHeadroomConverged), ...
                    'Kernel', 'on', 'Marker','.');
                set(h(1),'yscale','log');
                ylim([min(glutamateTimeToHeadroom(glutamateTimeToHeadroomConverged)) ...
                    max(glutamateTimeToHeadroom(glutamateTimeToHeadroomConverged))+max(glutamateTimeToHeadroom(glutamateTimeToHeadroomConverged))/5]);
                set(h(3),'xscale','log');
                title(['Glutamate Headroom vs Time to headroom',perturbationString(perturbation,0,1)]);
                xlabel('headroom'); ylabel('time to headroom');
                print([directory,glutamate_prep,'headroomVsTimeToHeadroom',perturbationString(perturbation,1,0)],'-dpng');
                print([directory,glutamate_prep,'headroomVsTimeToHeadroom',perturbationString(perturbation,1,0)],'-depsc');
                legend off;

                glutamateRiskyGroup = zeros(1,numel(glutamateHeadroomFinal));
                glutamateRiskyGroup(glutamateHeadroomFinal<headroomRisk ...
                    & glutamateTimeToHeadroom>=median(glutamateTimeToHeadroom)) = 1;
                glutamateRiskyGroup(glutamateHeadroomFinal>=headroomRisk ...
                    & glutamateTimeToHeadroom<median(glutamateTimeToHeadroom)) = 2; 
                glutamateRiskyGroup(glutamateHeadroomFinal>=headroomRisk ...
                    & glutamateTimeToHeadroom>=median(glutamateTimeToHeadroom)) = 3;

                [~, figNum] = newFigure(figNum, true);
                hp1 = uipanel('position', [0.0 0.5 0.5 0.5]);
                h = scatterhist(glutamateHeadroomFinal(glutamateTimeToHeadroomConverged), ...
                    glutamateTimeToHeadroom(glutamateTimeToHeadroomConverged), ...
                    'Group', glutamateRiskyGroup(glutamateTimeToHeadroomConverged), 'Kernel', 'on', 'Parent', hp1, ...
                    'Marker','.');
                set(h(1),'yscale','log');
                set(h(3),'xscale','log');
                ylim([min(glutamateTimeToHeadroom(glutamateTimeToHeadroomConverged)) ...
                    max(glutamateTimeToHeadroom(glutamateTimeToHeadroomConverged))+max(glutamateTimeToHeadroom(glutamateTimeToHeadroomConverged))/5]);
                xlabel('headroom'); ylabel('time to headroom');
    %                     legend off;

                hp2 = uipanel('position', [0.5 0.5 0.5 0.5]);
                tempHz = glutamateHzFinal(glutamatePreIndexs1D);
                scatterhist(spineAMPAWeights1D(glutamateTimeToHeadroomConverged), ...
                    tempHz(glutamateTimeToHeadroomConverged), ...
                    'Group', glutamateRiskyGroup(glutamateTimeToHeadroomConverged), 'Kernel', 'on', 'Parent', hp2, ...
                    'Marker','.');
                xlim([0 AMPAMax]); xlabel('glutamate bound to AMPA');
                ylim([0 HzMax]); 
                ylabel('in Hz');
                legend off;

                hp3 = uipanel('position', [0.0 0.0 0.5 0.5]);
                tempCleft = sum(glutamateCleftAstrocyteNeurotransmitter((measurementTrange(1)/sf)-(Trange(1)/sf): ...
                    (measurementTrange(end)/sf)-(Trange(1)/sf),:,:,:));
                scatterhist(spineAMPAWeights1D(glutamateTimeToHeadroomConverged), ...
                    tempCleft(glutamateTimeToHeadroomConverged), ...
                    'Group', glutamateRiskyGroup(glutamateTimeToHeadroomConverged), 'Kernel', 'on', 'Parent', hp3, ...
                    'Marker','.');
                xlim([0 AMPAMax]); xlabel('glutamate bound to AMPA');
                ylabel('glutamate flux');
                legend off;

                hp4 = uipanel('position', [0.5 0.0 0.5 0.5]);
                scatterhist(tempHz(glutamateTimeToHeadroomConverged), ...
                    tempCleft(glutamateTimeToHeadroomConverged), ...
                    'Group', glutamateRiskyGroup(glutamateTimeToHeadroomConverged), 'Kernel', 'on', 'Parent', hp4, ...
                    'Marker','.');
                xlim([0 HzMax]); xlabel('in Hz');
                ylabel('glutamate flux');
                legend off;

                print([directory,glutamate_prep,'riskySynapses',perturbationString(perturbation,1,0)],'-dpng');
                print([directory,glutamate_prep,'riskySynapses',perturbationString(perturbation,1,0)],'-depsc');
                clear tempHz;

                [~, figNum] = newFigure(figNum, true);
                bar([0,1,2,3],[mean(tempCleft(glutamateRiskyGroup==0 & glutamateTimeToHeadroomConverged)), ...
                    mean(tempCleft(glutamateRiskyGroup==1 & glutamateTimeToHeadroomConverged)), ...
                    mean(tempCleft(glutamateRiskyGroup==2 & glutamateTimeToHeadroomConverged)), ...
                    mean(tempCleft(glutamateRiskyGroup==3 & glutamateTimeToHeadroomConverged))]);
                title('Glutamate flux for each risky group');
                xlabel('risky group');
                ylabel('mean glutamate flux');

                print([directory,glutamate_prep,'riskySynapsesCleftGlutamate',perturbationString(perturbation,1,0)],'-dpng');
                print([directory,glutamate_prep,'riskySynapsesCleftGlutamate',perturbationString(perturbation,1,0)],'-depsc');
                clear tempCleft;
            end
            if (postprocess_GABA)
                [~, figNum] = newFigure(figNum, true);
                h = scatterhist(GABAHeadroomFinal(GABATimeToHeadroomConverged), ...
                    GABATimeToHeadroom(GABATimeToHeadroomConverged), ...
                    'Kernel', 'on', 'Marker','.');
                set(h(1),'yscale','log');
                ylim([min(GABATimeToHeadroom(GABATimeToHeadroomConverged)) ...
                    max(GABATimeToHeadroom(GABATimeToHeadroomConverged))+max(GABATimeToHeadroom(GABATimeToHeadroomConverged))/5]);
                set(h(3),'xscale','log');
                title(['GABA Headroom vs Time to headroom',perturbationString(perturbation,0,1)]);
                xlabel('headroom'); ylabel('time to headroom');
                print([directory,GABA_prep,'headroomVsTimeToHeadroom',perturbationString(perturbation,1,0)],'-dpng');
                print([directory,GABA_prep,'headroomVsTimeToHeadroom',perturbationString(perturbation,1,0)],'-depsc');
                legend off;

                GABARiskyGroup = zeros(1,numel(GABAHeadroomFinal));
                GABARiskyGroup(GABAHeadroomFinal<headroomRisk ...
                    & GABATimeToHeadroom>=median(GABATimeToHeadroom)) = 1;
                GABARiskyGroup(GABAHeadroomFinal>=headroomRisk ...
                    & GABATimeToHeadroom<median(GABATimeToHeadroom)) = 2; 
                GABARiskyGroup(GABAHeadroomFinal>=headroomRisk ...
                    & GABATimeToHeadroom>=median(GABATimeToHeadroom)) = 3;

                [~, figNum] = newFigure(figNum, true);
                hp1 = uipanel('position', [0.0 0.5 0.5 0.5]);
                h = scatterhist(GABAHeadroomFinal(GABATimeToHeadroomConverged), ...
                    GABATimeToHeadroom(GABATimeToHeadroomConverged), ...
                    'Group', GABARiskyGroup(GABATimeToHeadroomConverged), 'Kernel', 'on', 'Parent', hp1, ...
                    'Marker','.');
                set(h(1),'yscale','log');
                set(h(3),'xscale','log');
                ylim([min(GABATimeToHeadroom(GABATimeToHeadroomConverged)) ...
                    max(GABATimeToHeadroom(GABATimeToHeadroomConverged))+max(GABATimeToHeadroom(GABATimeToHeadroomConverged))/5]);
                xlabel('headroom'); ylabel('time to headroom');
    %                     legend off;

                hp2 = uipanel('position', [0.5 0.5 0.5 0.5]);
                tempHz = GABAHzFinal(GABAPreIndexs1D);
                scatterhist(spineAMPAWeights1D(GABATimeToHeadroomConverged), ...
                    tempHz(GABATimeToHeadroomConverged), ...
                    'Group', GABARiskyGroup(GABATimeToHeadroomConverged), 'Kernel', 'on', 'Parent', hp2, ...
                    'Marker','.');
                xlim([0 AMPAMax]); xlabel('GABA bound to AMPA');
                ylim([0 HzMax]); 
                ylabel('in Hz');
                legend off;

                hp3 = uipanel('position', [0.0 0.0 0.5 0.5]);
                tempCleft = sum(GABACleftAstrocyteNeurotransmitter((measurementTrange(1)/sf)-(Trange(1)/sf): ...
                    (measurementTrange(end)/sf)-(Trange(1)/sf),:,:,:));
                scatterhist(spineAMPAWeights1D(GABATimeToHeadroomConverged), ...
                    tempCleft(GABATimeToHeadroomConverged), ...
                    'Group', GABARiskyGroup(GABATimeToHeadroomConverged), 'Kernel', 'on', 'Parent', hp3, ...
                    'Marker','.');
                xlim([0 AMPAMax]); xlabel('GABA bound to AMPA');
                ylabel('GABA flux');
                legend off;

                hp4 = uipanel('position', [0.5 0.0 0.5 0.5]);
                scatterhist(tempHz(GABATimeToHeadroomConverged), ...
                    tempCleft(GABATimeToHeadroomConverged), ...
                    'Group', GABARiskyGroup(GABATimeToHeadroomConverged), 'Kernel', 'on', 'Parent', hp4, ...
                    'Marker','.');
                xlim([0 HzMax]); xlabel('in Hz');
                ylabel('GABA flux');
                legend off;

                print([directory,GABA_prep,'riskySynapses',perturbationString(perturbation,1,0)],'-dpng');
                print([directory,GABA_prep,'riskySynapses',perturbationString(perturbation,1,0)],'-depsc');
                clear tempHz;

                [~, figNum] = newFigure(figNum, true);
                bar([0,1,2,3],[mean(tempCleft(GABARiskyGroup==0 & GABATimeToHeadroomConverged)), ...
                    mean(tempCleft(GABARiskyGroup==1 & GABATimeToHeadroomConverged)), ...
                    mean(tempCleft(GABARiskyGroup==2 & GABATimeToHeadroomConverged)), ...
                    mean(tempCleft(GABARiskyGroup==3 & GABATimeToHeadroomConverged))]);
                title('GABA flux for each risky group');
                xlabel('risky group');
                ylabel('mean GABA flux');

                print([directory,GABA_prep,'riskySynapsesCleftGABA',perturbationString(perturbation,1,0)],'-dpng');
                print([directory,GABA_prep,'riskySynapsesCleftGABA',perturbationString(perturbation,1,0)],'-depsc');
                clear tempCleft;
            end
        end
        if (postprocess_EIratio)
            %%
            glutamateE = zeros(numel(1:T(end)/sf),1);
            GABAI = zeros(numel(1:T(end)/sf),1);
            for t=1:T(end)/sf
                glutamateE(t) = mean(bsxfun(@min, glutamateAvailableNeurotransmitter(t, :), ...
                    spineAMPAWeights1D'));
                GABAI(t) = mean(bsxfun(@min, GABAAvailableNeurotransmitter(t, :), ...
                    dendriteGABAWeights1D'));
            end
            [~, figNum] = newFigure(figNum, false);
            subplot(2,1,1);
            hold on;
            plot(glutamateE);
            plot(GABAI);
            subplot(2,1,2);
            plot(glutamateE ./ GABAI);
            
            print([directory,'EI'],'-dpng');
            print([directory,'EI'],'-depsc');
            %%
        end
    end
    %%
    if (postprocess_Perturbation)
        if (~perturbation)
            % not perturbation yet so store previous data needed
            if (postprocess_PerturbationHz)
                beforePerturbationHz = HzFinal;
            end
            if (postprocess_PerturbationAMPA)
                beforePerturbationAMPAWeights1D = spineAMPAWeights1D;
            end
            beforeHeadroomFinal = headroomFinal;
            beforeTimeToHeadroom = timeToHeadroom;
            beforeTimeToHeadroomConverged = timeToHeadroomConverged;
            beforeRiskyGroup = riskyGroup;
        else
            afterHeadroomFinal = headroomFinal;
            afterTimeToHeadroom = timeToHeadroom;
            afterTimeToHeadroomConverged = timeToHeadroomConverged;
            clear headroomFinal timeToHeadroom timeToHeadroomConverged;
            diffHeadroomFinal = afterHeadroomFinal - ...
                beforeHeadroomFinal;
            diffTimeToHeadroom = afterTimeToHeadroom - ...
                beforeTimeToHeadroom;
            headroomConverged = beforeTimeToHeadroomConverged & ...
                afterTimeToHeadroomConverged;
            
            % now perturbation, so process
            for p=1:3 % 1=Hz, 2=AMPA, 3=both
                if (p==1)
                    if (postprocess_PerturbationHz)
                        afterPerturbationHz = HzFinal;
                        diffPerturbationHz = afterPerturbationHz - ...
                            beforePerturbationHz;
                        perturbationType = 'Hz';
                        clear HzFinal;
    
                        beforePerturbation = beforePerturbationHz(preIndexs1D);
                        afterPerturbation = afterPerturbationHz(preIndexs1D);
                        diffPerturbation = diffPerturbationHz(preIndexs1D);
                    else
                        continue;
                    end
                elseif (p==2)
                    if (postprocess_PerturbationAMPA)
                        afterPerturbationAMPAWeights1D = spineAMPAWeights1D;
                        diffPerturbationAMPA = afterPerturbationAMPAWeights1D - ...
                            beforePerturbationAMPAWeights1D;
                        perturbationType = 'AMPA';
                        clear AMPAWeights1D;
    
                        beforePerturbation = beforePerturbationAMPAWeights1D;
                        afterPerturbation = afterPerturbationAMPAWeights1D;
                        diffPerturbation = diffPerturbationAMPA;
                    else
                        continue;
                    end
                elseif (p==3)
                    if (postprocess_PerturbationHz && postprocess_PerturbationAMPA)
                        perturbationType = 'Both';
                    else
                        continue;
                    end
                else
                    continue;
                end
                
                % Hz or AMPA perturbations separately
                if (p < 3)
                    %% Plot perturbation
                    [~, figNum] = newFigure(figNum, true);
                    subplot(2,2,1);
                    bar(diffPerturbation);
                    xlim([0 numel(diffPerturbation)]);
                    title([perturbationType,' perturbation']);
                    xlabel('pre-synaptic afferent'); ylabel('perturbation');
                    subplot(2,2,2);
                    histogram(diffPerturbation);
                    title([perturbationType,' perturbation']);
                    xlabel('perturbation'); ylabel('count');
                    subplot(2,2,3);
                    scatter(beforePerturbation, afterPerturbation,'.');
                    if (p==1)
                        xlim([0 HzMax]); ylim([0 HzMax]);
                    elseif (p==2)
                        xlim([0 AMPAMax]); ylim([0 AMPAMax]);
                    end
                    h=lsline;
                    set(h,'color','r');
                    R=corrcoef(beforePerturbation, afterPerturbation);
                    R_squared=R(2)^2;
                    title([perturbationType,...
                        ' before and after perturbation ', '(R squared= ', ...
                        num2str(R_squared), ')']);
                    xlabel(['before ',perturbationType]);
                    ylabel(['after ',perturbationType]);

                    print([directory,'perturbation',perturbationType],'-dpng');
                    print([directory,'perturbation',perturbationType],'-depsc');
                    clear h R R_squared;
                    %% Compare before and after
                    [~, figNum] = newFigure(figNum, true);

                    subplot(2,3,1);
                    scatter(beforeHeadroomFinal(headroomConverged), ...
                        afterHeadroomFinal(headroomConverged),'.');
                    h=lsline;
                    xlim([0 AMPAMax]); ylim([0 AMPAMax]);
                    set(h,'color','r');
                    R=corrcoef(beforeHeadroomFinal(headroomConverged), ...
                        afterHeadroomFinal(headroomConverged));
                    R_squared=R(2)^2;
                    title({['Headroom before and after ',perturbationType,...
                        ' perturbation'], ...
                        ['(R squared= ', num2str(R_squared), ')']});
                    xlabel('before'); ylabel('after');

                    subplot(2,3,4);
                    scatter(beforeTimeToHeadroom(headroomConverged), ...
                        afterTimeToHeadroom(headroomConverged),'.');
                    set(gca,'xscale','log');
                    set(gca,'yscale','log');
                    h=lsline;
                    xlim([min(beforeTimeToHeadroom(headroomConverged)) ...
                        max(beforeTimeToHeadroom(headroomConverged))]);
                    ylim([min(afterTimeToHeadroom(headroomConverged)) ...
                        max(afterTimeToHeadroom(headroomConverged))]);
                    set(h,'color','r');
                    R=corrcoef(beforeTimeToHeadroom(headroomConverged), ...
                        afterTimeToHeadroom(headroomConverged));
                    R_squared=R(2)^2;
                    title({['Time to headroom before and after ',perturbationType,...
                        'perturbation'], ...
                        ['(R squared= ', num2str(R_squared), ')']});
                    xlabel('before'); ylabel('after');

                    subplot(2,3,2);
                    scatter(diffPerturbation(headroomConverged), ...
                        afterHeadroomFinal(headroomConverged),'.'); 
                    h=lsline;
                    if (p==1)
                        xlim([-HzMax HzMax]); 
                    elseif(p==2)
                        xlim([-AMPAMax AMPAMax]); 
                    end
                    ylim([0 AMPAMax]);
                    set(h,'color','r');
                    R=corrcoef(diffPerturbation(headroomConverged), ...
                        afterHeadroomFinal(headroomConverged));
                    R_squared=R(2)^2;
                    title({[perturbationType,' perturbation vs headroom'], ...
                        ['(R squared= ', num2str(R_squared), ')']});
                    xlabel([perturbationType,' perturbation']); 
                    ylabel('after');

                    subplot(2,3,5);
                    scatter(diffPerturbation(headroomConverged), ...
                        afterTimeToHeadroom(headroomConverged),'.');
                    set(gca,'yscale','log');
                    h=lsline;
                    if (p==1)
                        xlim([-HzMax HzMax]);
                    elseif (p==2)
                        xlim([-AMPAMax AMPAMax]);
                    end
                    ylim([min(afterTimeToHeadroom(headroomConverged)) ...
                        max(afterTimeToHeadroom(headroomConverged))]);
                    set(h,'color','r');
                    R=corrcoef(diffPerturbation(headroomConverged), ...
                        afterTimeToHeadroom(headroomConverged));
                    R_squared=R(2)^2;
                    title({[perturbationType,' perturbation vs time to headroom'], ...
                        ['(R squared= ', num2str(R_squared), ')']});  
                    xlabel([perturbationType,' perturbation']);
                    ylabel('after');    

                    subplot(2,3,3);
                    scatter(diffPerturbation(headroomConverged), ...
                        diffHeadroomFinal(headroomConverged),'.'); 
                    h=lsline;
                    if (p==1)
                        xlim([-HzMax HzMax]); 
                    elseif (p==2)
                        xlim([-AMPAMax AMPAMax]); 
                    end
                    ylim([-AMPAMax AMPAMax]);
                    set(h,'color','r');
                    R=corrcoef(diffPerturbation(headroomConverged), ...
                        diffHeadroomFinal(headroomConverged));
                    R_squared=R(2)^2;
                    title({[perturbationType,' perturbation vs headroom difference '], ...
                        ['(R squared= ', num2str(R_squared), ')']});
                    xlabel([perturbationType,' perturbation']);
                    ylabel('headroom difference');

                    subplot(2,3,6);
                    scatter(diffPerturbation(headroomConverged), ...
                        diffTimeToHeadroom(headroomConverged),'.');
                    set(gca,'yscale','log');
                    h=lsline;
                    if (p==1)
                        xlim([-HzMax HzMax]);
                    elseif (p==2)
                        xlim([-AMPAMax AMPAMax]);
                    end
                    ylim([min(diffTimeToHeadroom(headroomConverged)) ...
                        max(diffTimeToHeadroom(headroomConverged))]);
                    set(h,'color','r');
                    R=corrcoef(diffPerturbation(headroomConverged), ...
                        diffTimeToHeadroom(headroomConverged));
                    R_squared=R(2)^2;
                    title({[perturbationType,' perturbation vs time to ',...
                        'headroom difference'], ...
                        ['(R squared= ', num2str(R_squared), ')']});  
                    xlabel([perturbationType,' perturbation']);
                    ylabel('time to headroom difference'); 

                    print([directory,'perturbation',perturbationType,'Comparison'],'-dpng');
                    print([directory,'perturbation',perturbationType,'Comparison'],'-depsc');
                    clear h R R_squared;
                    %% Glutamate flux as a function of perturbation
                    [~, figNum] = newFigure(figNum, true);
                    
                    [glutamateFluxPerturbation, glutamateFluxPerturbationMean, ...
                        glutamateFluxPerturbationMeanBinned, ...
                        glutamateFluxPerturbationP] = ...
                        glutamateFluxVsPerturbation(cleftGlutamate, ...
                        afterTimeToHeadroom, diffPerturbation, ...
                        headroomConverged, p, HzMax, AMPAMax, ...
                        perturbationType, '');

                    print([directory,'perturbation',perturbationType,...
                        'GlutamateFlux'],'-dpng');
                    print([directory,'perturbation',perturbationType,...
                        'GlutamateFlux'],'-depsc');
                    
                    glutamateFluxPerturbationRisk = cell(1,4);
                    glutamateFluxPerturbationMeanRisk = cell(1,4);
                    glutamateFluxPerturbationMeanBinnedRisk = cell(1,4);
                    glutamateFluxPerturbationPRisk = cell(1,4);
                    for group=0:3
                        [~, figNum] = newFigure(figNum, true);
                    
                        [glutamateFluxPerturbationRisk{group+1}, ...
                            glutamateFluxPerturbationMeanRisk{group+1}, ...
                            glutamateFluxPerturbationMeanBinnedRisk{group+1}, ...
                            glutamateFluxPerturbationPRisk{group+1}] = ...
                            glutamateFluxVsPerturbation(...
                            cleftGlutamate(:,beforeRiskyGroup==group), ...
                            afterTimeToHeadroom(beforeRiskyGroup==group), ...
                            diffPerturbation(beforeRiskyGroup==group), ...
                            headroomConverged(beforeRiskyGroup==group), ...
                            p, HzMax, AMPAMax, perturbationType, ...
                            [' (risk group ', num2str(group-1),')']);

                        print([directory,'perturbation',perturbationType,...
                            'GlutamateFlux_riskGroup',num2str(group)],'-dpng');
                        print([directory,'perturbation',perturbationType,...
                            'GlutamateFlux_riskGroup',num2str(group)],'-depsc');
                    end
                end
                % Hz and AMPA perturbations together
                if (p==3)
                    % Correlation of Hz and AMPA perturbations
                    [~, figNum] = newFigure(figNum, false);
                    scatter(diffPerturbationHz(preIndexs1D), ...
                        diffPerturbationAMPA, '.');
                    xlim([-HzMax HzMax]); ylim([-AMPAMax AMPAMax]);
                    title('Hz perturbation vs AMPA perturbation');
                    xlabel('Hz perturbation');
                    ylabel('AMPA perturbation');
                    h=lsline;
                    set(h,'color','r');
                    R=corrcoef(diffPerturbationHz(preIndexs1D), ...
                        diffPerturbationAMPA);
                    R_squared=R(2)^2;
                    title([perturbationType,...
                        ' difference ', '(R squared= ', ...
                        num2str(R_squared), ')']);                    
                    print([directory,'perturbation',perturbationType],'-dpng'); 
                    print([directory,'perturbation',perturbationType],'-depsc');                    
                    clear h R R_squared;
                    %%                    
                    HzRange = -140:20:140;
                    HzOffset = 10;
                    AMPARange = -1.4:0.2:1.4;
                    AMPAOffset = 0.1;
                    diffHeadroomBothPerturbationBinned = zeros(numel(HzRange),...
                        numel(AMPARange));
                    cleftGlutamateBothMeanPerturbationBinned = ...
                        zeros(numel(HzRange), numel(AMPARange));
                    cleftGlutamateBothMaxPerturbationBinned = ...
                        zeros(numel(HzRange), numel(AMPARange));
                    cleftGlutamateBothMeanPerturbationBinnedRisky = cell(1,4);
                    cleftGlutamateBothMaxPerturbationBinnedRisky = cell(1,4);
                    for i=1:4
                        cleftGlutamateBothMeanPerturbationBinnedRisky{i} = ...
                            zeros(numel(HzRange), numel(AMPARange));
                        cleftGlutamateBothMaxPerturbationBinnedRisky{i} = ...
                            zeros(numel(HzRange), numel(AMPARange));
                    end
                    clear i;
                    
                    for i=1:numel(HzRange)
                        for j=1:numel(AMPARange)
                            tempHz = diffPerturbationHz(preIndexs1D)...
                                    >=HzRange(i)-HzOffset & ...
                                diffPerturbationHz(preIndexs1D)...
                                    <HzRange(i)+HzOffset;
                            tempAMPA = diffPerturbationAMPA...
                                    >=AMPARange(j)-AMPAOffset & ...
                                diffPerturbationAMPA...
                                    <AMPARange(j)+AMPAOffset;
                            tempAll = tempHz & tempAMPA & headroomConverged';
                            diffHeadroomBothPerturbationBinned(i,j) = ...
                                mean(diffHeadroomFinal(tempAll));
                            cleftGlutamateBothMeanPerturbationBinned(i,j) = ...
                                mean(glutamateFluxPerturbation(tempAll));
                            temp = max(glutamateFluxPerturbation(tempAll));
                            if numel((temp) > 0)
                                cleftGlutamateBothMaxPerturbationBinned(i,j) = ...
                                    temp;
                            end
                            for group=0:3
                                tempAllRisky = tempAll & (beforeRiskyGroup==group)';
                                cleftGlutamateBothMeanPerturbationBinnedRisky{group+1}(i,j) = ...
                                    mean(glutamateFluxPerturbation(tempAllRisky));   
                                temp = max(glutamateFluxPerturbation(tempAllRisky));
                                if (numel(temp) > 0)
                                    cleftGlutamateBothMaxPerturbationBinnedRisky{group+1}(i,j) = ...
                                        temp;
                                end
                            end
                            clear group tempAllRisky;
                        end
                    end
                    
                    [~, figNum] = newFigure(figNum, false);
                    
                    subplot(2,2,1);
                    imagesc(diffHeadroomBothPerturbationBinned);
                    colormap(flipud(hot)); colorbar();
                    set(gca,'ydir','normal');
                    title('Headroom difference vs both perturbations');
                    xlabel('AMPA perturbation');
                    ylabel('Hz perturbation');
                    h = colorbar; ylabel(h, 'headroom difference')
                    yticks([1,3,5,7,8,9,11,13,15]);
                    xticks([1,3,5,7,8,9,11,13,15]);
                    xticklabels({'-1.4','-1.0','-0.6','-0.2','0.0',...
                        '0.2','0.6','0.8','1.0','1.4'});
                    yticklabels({'-140','-100','-60','-20','0',...
                        '20','60','100','140'});
                    
                    subplot(2,2,2);
                    imagesc(cleftGlutamateBothMeanPerturbationBinned);
                    colormap(flipud(hot)); colorbar();
                    set(gca,'ydir','normal');
                    title('Glutamate flux vs both perturbations');
                    xlabel('AMPA perturbation');
                    ylabel('Hz perturbation');
                    h = colorbar; ylabel(h, 'glutamate flux')
                    xticks([1,3,5,7,8,9,11,13,15]);
                    xticklabels({'-1.4','-1.0','-0.6','-0.2','0.0',...
                        '0.2','0.6','0.8','1.0','1.4'});   
                    yticks([1,3,5,7,8,9,11,13,15]);
                    yticklabels({'-140','-100','-60','-20','0',...
                        '20','60','100','140'});      
                    
                    subplot(2,2,3);
                    imagesc(cleftGlutamateBothMaxPerturbationBinned);
                    colormap(flipud(hot)); colorbar();
                    set(gca,'ydir','normal');
                    title('Glutamate peak vs both perturbations');
                    xlabel('AMPA perturbation');
                    ylabel('Hz perturbation');
                    h = colorbar; ylabel(h, 'glutamate peak')
                    xticks([1,3,5,7,8,9,11,13,15]);
                    xticklabels({'-1.4','-1.0','-0.6','-0.2','0.0',...
                        '0.2','0.6','0.8','1.0','1.4'});   
                    yticks([1,3,5,7,8,9,11,13,15]);
                    yticklabels({'-140','-100','-60','-20','0',...
                        '20','60','100','140'});               
                    
                    print([directory,'perturbation',perturbationType,'Comparison'],'-dpng');
                    print([directory,'perturbation',perturbationType,'Comparison'],'-depsc');
                    
                    [~, figNum] = newFigure(figNum, true);
                    
                    for group=1:4
                        subplot(2,2,group);
                        imagesc(cleftGlutamateBothMeanPerturbationBinnedRisky{group});
                        colormap(flipud(hot)); colorbar();
                        set(gca,'ydir','normal');
                        title(['Glutamate flux vs both perturbations (risk group ',...
                            num2str(group-1),')']);
                        xlabel('AMPA perturbation');
                        ylabel('Hz perturbation');
                        h = colorbar; ylabel(h, 'glutamate flux')
                        xticks([1,3,5,7,8,9,11,13,15]);
                        xticklabels({'-1.4','-1.0','-0.6','-0.2','0.0',...
                            '0.2','0.6','0.8','1.0','1.4'});   
                        yticks([1,3,5,7,8,9,11,13,15]);
                        yticklabels({'-140','-100','-60','-20','0',...
                            '20','60','100','140'});
                    end
                    
                    print([directory,'perturbation',perturbationType,'ComparisonRiskyMean'],'-dpng');
                    print([directory,'perturbation',perturbationType,'ComparisonRiskyMean'],'-depsc');
                    
                    [~, figNum] = newFigure(figNum, true);
                    
                    for group=1:4
                        subplot(2,2,group);
                        imagesc(cleftGlutamateBothMaxPerturbationBinnedRisky{group});
                        colormap(flipud(hot)); colorbar();
                        set(gca,'ydir','normal');
                        title(['Glutamate peak vs both perturbations (risk group ',...
                            num2str(group-1),')']);
                        xlabel('AMPA perturbation');
                        ylabel('Hz perturbation');
                        h = colorbar; ylabel(h, 'glutamate peak')
                        xticks([1,3,5,7,8,9,11,13,15]);
                        xticklabels({'-1.4','-1.0','-0.6','-0.2','0.0',...
                            '0.2','0.6','0.8','1.0','1.4'});   
                        yticks([1,3,5,7,8,9,11,13,15]);
                        yticklabels({'-140','-100','-60','-20','0',...
                            '20','60','100','140'});
                    end
                    
                    print([directory,'perturbation',perturbationType,'ComparisonRiskyMax'],'-dpng');
                    print([directory,'perturbation',perturbationType,'ComparisonRiskyMax'],'-depsc');
                    
                    clear HzRange HzOffset AMPARange AMPAOffset i j ...
                        tempHz tempAMPA tempCleft tempAll;
                end
            end
        end
    end
end
%% Functions
function [h, figNum] = newFigure(figNum, maximize)
    figNum = figNum + 1;
    h = figure(figNum);
    clf;
    if (maximize)
        set(h, 'Position', get(0,'Screensize'));
    end
end
function [ret] = perturbationString(perturbation, file, label)
    if (perturbation)
        if (file)
            ret = '_perturbation';
        elseif (label)
            ret = ' (Perturbation) ';
        end
    else
        ret = '';
    end
end
function [ret] = loadSpikes(directory, file, fileExt, ...
    spikesFilter, spikesFilterTmin, spikesFilterTmax, Tmin, Tmax, dt, ...
    XdimIn, YdimIn, ZdimIn)
    fid = fopen([directory,file,fileExt],'r');
    Xdim = fread(fid, 1, 'int');
    Ydim = fread(fid, 1, 'int');
    Zdim = fread(fid, 1, 'int');
    % {3D cell: X, Y, Z} (1D: spike times)
    ret = cell(Xdim,Ydim,Zdim);
    if ((Xdim ~= XdimIn) || (Ydim ~= YdimIn) || (Zdim ~= ZdimIn))
        disp('Dimensions wrong');
        return;
    end
    % temp(1,:) is neuron id and temp(2,:) is spike time
    temp = fread(fid, Inf, 'int'); % have to load all sadly
    if (mod(numel(temp),2) ~= 0)
        temp(end) = [];
    end
    temp = reshape(temp, [2, numel(temp)/2]);
    fclose(fid);
    clear fid;
    if (spikesFilter)
        temp(:,temp(2,:)<spikesFilterTmin/dt) = []; % filter out less than Tmin
        temp(:,temp(2,:)>spikesFilterTmax/dt) = []; % filter out greater than Tmax
    else
        temp(:,temp(2,:)<Tmin/dt) = []; % filter out greater than Tmin
        temp(:,temp(2,:)>Tmax/dt) = []; % filter out greater than Tmax
    end
    temp(1,:) = temp(1,:)+1; % because of 0 indexing
    temp = sortrows(temp'); % Soft by neuron id
    [~, temp_i] = unique(temp(:,1)); % First spike of each neuron
    spikesN = circshift(temp_i,-1) - temp_i; % Number of spikes for each neuron
    spikesN(end) = size(temp,1) - sum(spikesN(1:end-1)); % fix last neuron
    i=1; % Neuron id counter
    k=1; % Position in temp_i and spikesN
    for z=1:Zdim
        for y=1:Ydim
            for x=1:Xdim
                % If this neuron has spikes
                if (temp(temp_i(k),1) == i)
                    ret{x,y,z} = ...
                        temp(temp_i(k):temp_i(k)+spikesN(k)-1,2);
                    k=k+1;
                end
                i=i+1;
                if (k > size(unique(temp(:,1)),1))
                    break;
                end
            end
            if (k > size(unique(temp(:,1)),1))
                break;
            end
        end
        if (k > size(unique(temp(:,1)),1))
            break;
        end
    end
end               
function [ret] = load4D(directory, file, ...
    fileExt, Tmin, Tmax, sf, ...
    XdimInnerIn, YdimInnerIn, ZdimInnerIn)
    % Load from file
    fid = fopen([directory,file,fileExt],'r');
    XdimInner = fread(fid, 1, 'int');
    YdimInner = fread(fid, 1, 'int');
    ZdimInner = fread(fid, 1, 'int');
    dimAdjust = false;
    if ((XdimInner ~= XdimInnerIn) || (YdimInner ~= YdimInnerIn) ...
            || (ZdimInner ~= ZdimInnerIn))
        disp('Dimensions different ...');
        dimAdjust = true;
    end
    fseek(fid, ((XdimInner*YdimInner*ZdimInner)*(Tmin/sf))*4, 'cof');
    ret = fread(fid, (XdimInner*YdimInner*ZdimInner)*((Tmax-Tmin)/sf), 'float');
    fclose(fid);
    clear fid;
    % 4D: time, X, Y, Z
    ret = reshape(ret(1:floor(numel(ret)/...
        ((XdimInner*YdimInner*ZdimInner)))*XdimInner*YdimInner*ZdimInner), ...
        [XdimInner, YdimInner, ZdimInner, ...
        floor(numel(ret)/(XdimInner*YdimInner*ZdimInner))]);
    ret = permute(ret, [4, 1, 2, 3]);
    % Adjust if needed due to dimension differences
    if (dimAdjust)
        if (XdimInnerIn < XdimInner)
            disp('... adjusting X dimension');
            ret(:,XdimInnerIn+1:end,:,:) = [];
        end
        if (YdimInnerIn < YdimInner)
            disp('... adjusting Y dimension');
            ret(:,:,YdimInnerIn+1:end,:) = [];
        end
        if (ZdimInnerIn < ZdimInner)
            disp('... adjusting Z dimension');
            ret(:,:,:,ZdimInnerIn+1:end) = [];
        end
    end
end
function [ret, ret1D] = load2D(directory, file, fileExt, XdimInner, ...
    YdimInner, ZdimInner)
    fid = fopen([directory,file,fileExt],'r');
    % 4D: [1=pre 2=post], X, Y, Z
    ret = fread(fid, [2, Inf], 'int');
    ret1D = ret(1,:,:,:);
    fclose(fid);
    clear fid;
    ret = reshape(ret, [2, XdimInner, YdimInner, ZdimInner]);
end
function [ret, ret1D] = load3D(directory, file, fileExt, XdimInner, ...
    YdimInner, ZdimInner)
    fid = fopen([directory,file,fileExt],'r');
    % 3D: X, Y, Z
    ret1D = fread(fid, Inf, 'float');
    fclose(fid);
    clear fid;
    ret = reshape(ret1D, [XdimInner, YdimInner, ZdimInner]);
end
function plotModulation(directory, file, fileExt, titleStr, x, ...
    y, saveFile, fontsize)
    % mGluR5 modulation function only first
    fid = fopen([directory,file,fileExt],'r');
    mGluR5modulation = fread(fid, Inf, 'float');
    fclose(fid);
    clear fid;
    plot(0:1/1000:2, mGluR5modulation);
    title(titleStr);
    xlabel(x); ylabel(y);
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    print([directory,saveFile],'-dpng');
    print([directory,saveFile],'-depsc');
end
function [ret] = calculateHz(input, Tmin, Tmax, HzSf, Xdim, ...
    Ydim, Zdim, dt)
    ret=zeros(numel(Tmin:HzSf:Tmax-HzSf),Xdim,Ydim,Zdim);
    for x=1:Xdim
        for y=1:Ydim
            for z=1:Zdim
                i=1;
                for t=Tmin:HzSf:Tmax-HzSf
                   ret(i,x,y,z) = ...
                        numel(find(...
                            input{x,y,z}>=t/dt & ...
                            input{x,y,z}<(t+1)/dt...
                        ))*(1.0/HzSf);
                    i=i+1;                                
                end
            end
        end
    end
end
function [ret] = hist3D(range, in)
    ret = zeros(numel(range),size(in,1));
    for i=1:size(in,1)
        temp = in(i,:,:,:);
        ret(:,i) = hist(temp(:),range);
    end
end
function [ret] = histHeadroom3D(range, availableGlutamate, weights)
    ret = zeros(numel(range),size(availableGlutamate,1));  
    for i=1:size(availableGlutamate,1) % over time, glutamate and AMPA are same size
        % calculate the diff/headroom for each synapse and hist it
        tempAvailableGlutamate = availableGlutamate(i,:,:,:);
        ret(:,i) = hist(tempAvailableGlutamate(:)-weights, ...
            range);
    end
end
function [...
    glutamateHeadroomFinal, glutamateHeadroomFinalStd, ...
    glutamateHzFinal, glutamateTimeToHeadroom, glutamateTimeToHeadroomConverged, ...
    glutamateCB1RFinal, glutamateAvailableNeurotransmitterFinal, spineCaFinal, ...
    glutamateCleftAstrocyteeCBFinal, ...
    ...
    GABAHeadroomFinal, GABAHeadroomFinalStd, ...
    GABAHzFinal, GABATimeToHeadroom, GABATimeToHeadroomConverged, ...
    GABACB1RFinal, GABAAvailableNeurotransmitterFinal, dendriteCaFinal, ...
    GABACleftAstrocyteeCBFinal, ...
    ...
    outHzFinal, ...
    ...
    figNum] = ...
    plotComponentVsHeadroom(figNum, Xdim, Ydim, Zdim, ...
    XdimInner, YdimInner, ZdimInner, ...
    ...
    glutamateAvailableNeurotransmitter, AMPAWeights, AMPAWeights1D, glutamateInHz, ...
    glutamatePreIndexs1D, AMPAMax, glutamateCB1R, spineCa, ...
    glutamateCleftAstrocyteeCB, OutputAMPAIndexs1D, ...
    ...
    GABAAvailableNeurotransmitter, GABAWeights, GABAWeights1D, GABAInHz, ...
    GABAPreIndexs1D, GABAMax, GABACB1R, dendriteCa, ...
    GABACleftAstrocyteeCB, OutputGABARIndexs1D, ...
    ...
    outHz, ...
    ...
    finalHeadroomT, finalHzT, finalHeadroomStd, Tmin, sf, HzSf, ...
    inputSpikeTimesTmin, HzMax, ...
    directory, glutamate_prep, GABA_prep, fontsize)

    %% Calculate glutamate final headroom and time to headroom    
    glutamateHeadroomFinal = zeros(1,XdimInner*YdimInner*ZdimInner);
    glutamateHeadroomFinalStd = zeros(1,XdimInner*YdimInner*ZdimInner);
    glutamateTimeToHeadroom = zeros(1,XdimInner*YdimInner*ZdimInner);
    i = 1;
    for x=1:XdimInner
        for y=1:YdimInner
            for z=1:ZdimInner
                tempAvailableGlutamate = glutamateAvailableNeurotransmitter(:,x,y,z);
                % Headroom
                glutamateHeadroomFinal(i) = mean(tempAvailableGlutamate( ...
                    ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf) - ...
                    AMPAWeights(x,y,z));
                glutamateHeadroomFinalStd(i) = std(tempAvailableGlutamate( ...
                    ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf) - ...
                    AMPAWeights(x,y,z));
                % Time to headroom
                tempHeadroom = tempAvailableGlutamate - AMPAWeights(x,y,z);
                temp = find(tempHeadroom<=glutamateHeadroomFinal(i)...
                    +glutamateHeadroomFinalStd(i)*(finalHeadroomStd/2) & ...
                    tempHeadroom>=glutamateHeadroomFinal(i)...
                    -glutamateHeadroomFinalStd(i)*(finalHeadroomStd/2));
                if (numel(temp) == 0)
                    glutamateTimeToHeadroom(i) = -9;
                else
                    glutamateTimeToHeadroom(i) = temp(1);
                end
                i = i + 1;
            end
        end
    end
    glutamateTimeToHeadroomConverged = glutamateTimeToHeadroom>1; % those that actually had to converge
    clear tempAvailableGlutamate temp;
    
    %% Calculate GABA final headroom and time to headroom    
    GABAHeadroomFinal = zeros(1,XdimInner*YdimInner*ZdimInner);
    GABAHeadroomFinalStd = zeros(1,XdimInner*YdimInner*ZdimInner);
    GABATimeToHeadroom = zeros(1,XdimInner*YdimInner*ZdimInner);
    i = 1;
    for x=1:XdimInner
        for y=1:YdimInner
            for z=1:ZdimInner
                tempAvailableGABA = GABAAvailableNeurotransmitter(:,x,y,z);
                % Headroom
                GABAHeadroomFinal(i) = mean(tempAvailableGABA( ...
                    ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf) - ...
                    GABAWeights(x,y,z));
                GABAHeadroomFinalStd(i) = std(tempAvailableGABA( ...
                    ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf) - ...
                    GABAWeights(x,y,z));
                % Time to headroom
                tempHeadroom = tempAvailableGABA - GABAWeights(x,y,z);
                temp = find(tempHeadroom<=GABAHeadroomFinal(i)...
                    +GABAHeadroomFinalStd(i)*(finalHeadroomStd/2) & ...
                    tempHeadroom>=GABAHeadroomFinal(i)...
                    -GABAHeadroomFinalStd(i)*(finalHeadroomStd/2));
                if (numel(temp) == 0)
                    GABATimeToHeadroom(i) = -9;
                else
                    GABATimeToHeadroom(i) = temp(1);
                end
                i = i + 1;
            end
        end
    end
    GABATimeToHeadroomConverged = GABATimeToHeadroom>1; % those that actually had to converge
    clear tempAvailableGABA temp i x y z;
    
    %% Calculate mean final values
    glutamateHzFinal = zeros(Xdim*Ydim*Zdim,1);
    GABAHzFinal = zeros(Xdim*Ydim*Zdim,1);
    outHzFinal = zeros(Xdim*Ydim*Zdim,1);
    i=1;
    for x=1:Xdim
        for y=1:Ydim
            for z=1:Zdim
                glutamateHzFinal(i) = mean(glutamateInHz(((finalHzT(1)-inputSpikeTimesTmin)/HzSf)+1: ...
                    (finalHzT(end)-inputSpikeTimesTmin)/HzSf,x,y,z));
                GABAHzFinal(i) = mean(GABAInHz(((finalHzT(1)-inputSpikeTimesTmin)/HzSf)+1: ...
                    (finalHzT(end)-inputSpikeTimesTmin)/HzSf,x,y,z));
                outHzFinal(i) = mean(outHz(((finalHzT(1)-inputSpikeTimesTmin)/HzSf)+1: ...
                    (finalHzT(end)-inputSpikeTimesTmin)/HzSf,x,y,z));
                i=i+1;
            end
        end
    end 
    glutamateCB1RFinal = zeros(XdimInner*YdimInner*ZdimInner,1);
    GABACB1RFinal = zeros(XdimInner*YdimInner*ZdimInner,1);
    glutamateAvailableNeurotransmitterFinal = zeros(XdimInner*YdimInner*ZdimInner,1);
    GABAAvailableNeurotransmitterFinal = zeros(XdimInner*YdimInner*ZdimInner,1);
    spineCaFinal = zeros(XdimInner*YdimInner*ZdimInner,1);
    dendriteCaFinal = zeros(XdimInner*YdimInner*ZdimInner,1);
    glutamateCleftAstrocyteeCBFinal = zeros(XdimInner*YdimInner*ZdimInner,1);
    GABACleftAstrocyteeCBFinal = zeros(XdimInner*YdimInner*ZdimInner,1);
    i=1;
    for x=1:XdimInner
        for y=1:YdimInner
            for z=1:ZdimInner
                glutamateCB1RFinal(i) = mean(glutamateCB1R( ...
                    ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf,x,y,z));
                GABACB1RFinal(i) = mean(GABACB1R( ...
                    ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf,x,y,z));
                glutamateAvailableNeurotransmitterFinal(i) = mean(glutamateAvailableNeurotransmitter( ...
                    ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf,x,y,z));
                GABAAvailableNeurotransmitterFinal(i) = mean(GABAAvailableNeurotransmitter( ...
                    ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf,x,y,z));                
                spineCaFinal(i) = mean(spineCa( ...
                    ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf,x,y,z));
                dendriteCaFinal(i) = mean(dendriteCa( ...
                    ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf,x,y,z));                               
                glutamateCleftAstrocyteeCBFinal(i) = mean(glutamateCleftAstrocyteeCB( ...
                    ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf,x,y,z));
                GABACleftAstrocyteeCBFinal(i) = mean(GABACleftAstrocyteeCB( ...
                    ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf,x,y,z));              
%                 spineCaFinal(i) = sum(spineCa( ...
%                     ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf,x,y,z));
%                 dendriteCaFinal(i) = sum(dendriteCa( ...
%                     ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf,x,y,z));                               
%                 glutamateCleftAstrocyteeCBFinal(i) = sum(glutamateCleftAstrocyteeCB( ...
%                     ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf,x,y,z));
%                 GABACleftAstrocyteeCBFinal(i) = sum(GABACleftAstrocyteeCB( ...
%                     ((finalHeadroomT(1)-Tmin)/sf)+1:(finalHeadroomT(end)-Tmin)/sf,x,y,z));
                i = i + 1;
            end
        end
    end
    
    %% Figure properties
    [~, figNum] = newFigure(figNum, true);   
    glutamateHeadroomFigNum = figNum;
    [~, figNum] = newFigure(figNum, true);   
    glutamateTimeToHeadroomFigNum = figNum;
    [~, figNum] = newFigure(figNum, true);   
    GABAHeadroomFigNum = figNum;
    [~, figNum] = newFigure(figNum, true);   
    GABATimeToHeadroomFigNum = figNum;
    
    %% Glutamate headroom
    figure(glutamateHeadroomFigNum);
    
    subplot(2,7,1);    
    scatter(AMPAWeights1D, glutamateHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(AMPAWeights1D, glutamateHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('AMPA'); ylabel('Glutamate headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    subplot(2,7,8);    
    scatter(GABAWeights1D, glutamateHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABAWeights1D, glutamateHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA'); ylabel('Glutamate headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,2);
    scatter(glutamateHzFinal(glutamatePreIndexs1D), glutamateHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateHzFinal(glutamatePreIndexs1D), glutamateHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Glutamate Hz'); ylabel('Glutamate headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,9);
    scatter(GABAHzFinal(GABAPreIndexs1D), glutamateHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABAHzFinal(GABAPreIndexs1D), glutamateHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA Hz'); ylabel('Glutamate headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,3);
    scatter(glutamateCB1RFinal, glutamateHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateCB1RFinal, glutamateHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Glutamate CB1R'); ylabel('Glutamate headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,10);
    scatter(GABACB1RFinal, glutamateHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABACB1RFinal, glutamateHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA CB1R'); ylabel('Glutamate headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,4);
    scatter(glutamateAvailableNeurotransmitterFinal, glutamateHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateAvailableNeurotransmitterFinal, glutamateHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Available glutamate'); ylabel('Glutamate headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,11);
    scatter(GABAAvailableNeurotransmitterFinal, glutamateHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABAAvailableNeurotransmitterFinal, glutamateHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Available GABA'); ylabel('Glutamate headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,5);
    scatter(spineCaFinal, glutamateHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(spineCaFinal, glutamateHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Spine Ca2+'); ylabel('Glutamate headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,12);
    scatter(dendriteCaFinal, glutamateHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(dendriteCaFinal, glutamateHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Dendrite Ca2+'); ylabel('Glutamate headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,6);
    scatter(glutamateCleftAstrocyteeCBFinal, glutamateHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateCleftAstrocyteeCBFinal, glutamateHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Glutamate cleft eCB'); ylabel('Glutamate headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,13);
    scatter(GABACleftAstrocyteeCBFinal, glutamateHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABACleftAstrocyteeCBFinal, glutamateHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA cleft eCB'); ylabel('Glutamate headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,[7,14]);
    scatter(outHzFinal(OutputAMPAIndexs1D), glutamateHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(outHzFinal(OutputAMPAIndexs1D), glutamateHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Out Hz'); ylabel('Glutamate headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    print([directory,glutamate_prep,'components_vs_headroom'],'-dpng');
    print([directory,glutamate_prep,'components_vs_headroom'],'-depsc');
    
    %% Glutamate time to headroom
    figure(glutamateTimeToHeadroomFigNum);
    
    subplot(2,7,1);    
    scatter(AMPAWeights1D, glutamateTimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(AMPAWeights1D, glutamateTimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('AMPA'); ylabel('Glutamate time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    subplot(2,7,8);    
    scatter(GABAWeights1D, glutamateTimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABAWeights1D, glutamateTimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA'); ylabel('Glutamate time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,2);
    scatter(glutamateHzFinal(glutamatePreIndexs1D), glutamateTimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateHzFinal(glutamatePreIndexs1D), glutamateTimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Glutamate Hz'); ylabel('Glutamate time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,9);
    scatter(GABAHzFinal(GABAPreIndexs1D), glutamateTimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABAHzFinal(GABAPreIndexs1D), glutamateTimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA Hz'); ylabel('Glutamate time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,3);
    scatter(glutamateCB1RFinal, glutamateTimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateCB1RFinal, glutamateTimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Glutamate CB1R'); ylabel('Glutamate time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,10);
    scatter(GABACB1RFinal, glutamateTimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABACB1RFinal, glutamateTimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA CB1R'); ylabel('Glutamate time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,4);
    scatter(glutamateAvailableNeurotransmitterFinal, glutamateTimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateAvailableNeurotransmitterFinal, glutamateTimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Available glutamate'); ylabel('Glutamate time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,11);
    scatter(GABAAvailableNeurotransmitterFinal, glutamateTimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABAAvailableNeurotransmitterFinal, glutamateTimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Available GABA'); ylabel('Glutamate time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,5);
    scatter(spineCaFinal, glutamateTimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(spineCaFinal, glutamateTimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Spine Ca2+'); ylabel('Glutamate time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,12);
    scatter(dendriteCaFinal, glutamateTimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(dendriteCaFinal, glutamateTimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Dendrite Ca2+'); ylabel('Glutamate time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,6);
    scatter(glutamateCleftAstrocyteeCBFinal, glutamateTimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateCleftAstrocyteeCBFinal, glutamateTimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Glutamate cleft eCB'); ylabel('Glutamate time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,13);
    scatter(GABACleftAstrocyteeCBFinal, glutamateTimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABACleftAstrocyteeCBFinal, glutamateTimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA cleft eCB'); ylabel('Glutamate time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,[7,14]);
    scatter(outHzFinal(OutputAMPAIndexs1D), glutamateTimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(outHzFinal(OutputAMPAIndexs1D), glutamateTimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Out Hz'); ylabel('Glutamate time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    print([directory,glutamate_prep,'components_vs_timeToHeadroom'],'-dpng');
    print([directory,glutamate_prep,'components_vs_timeToHeadroom'],'-depsc');
    
    %% GABA headroom
    figure(GABAHeadroomFigNum);
    
    subplot(2,7,1);    
    scatter(AMPAWeights1D, GABAHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(AMPAWeights1D, GABAHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('AMPA'); ylabel('GABA headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    subplot(2,7,8);    
    scatter(GABAWeights1D, GABAHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABAWeights1D, GABAHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA'); ylabel('GABA headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,2);
    scatter(glutamateHzFinal(glutamatePreIndexs1D), GABAHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateHzFinal(glutamatePreIndexs1D), GABAHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Glutamate Hz'); ylabel('GABA headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,9);
    scatter(GABAHzFinal(GABAPreIndexs1D), GABAHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABAHzFinal(GABAPreIndexs1D), GABAHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA Hz'); ylabel('GABA headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,3);
    scatter(glutamateCB1RFinal, GABAHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateCB1RFinal, GABAHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Glutamate CB1R'); ylabel('GABA headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,10);
    scatter(GABACB1RFinal, GABAHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABACB1RFinal, GABAHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA CB1R'); ylabel('GABA headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,4);
    scatter(glutamateAvailableNeurotransmitterFinal, GABAHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateAvailableNeurotransmitterFinal, GABAHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Available glutamate'); ylabel('GABA headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,11);
    scatter(GABAAvailableNeurotransmitterFinal, GABAHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABAAvailableNeurotransmitterFinal, GABAHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Available GABA'); ylabel('GABA headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,5);
    scatter(spineCaFinal, GABAHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(spineCaFinal, GABAHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Spine Ca2+'); ylabel('GABA headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,12);
    scatter(dendriteCaFinal, GABAHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(dendriteCaFinal, GABAHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Dendrite Ca2+'); ylabel('GABA headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,6);
    scatter(glutamateCleftAstrocyteeCBFinal, GABAHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateCleftAstrocyteeCBFinal, GABAHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Glutamate cleft eCB'); ylabel('GABA headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,13);
    scatter(GABACleftAstrocyteeCBFinal, GABAHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABACleftAstrocyteeCBFinal, GABAHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA cleft eCB'); ylabel('GABA headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,[7,14]);
    scatter(outHzFinal(OutputGABARIndexs1D), GABAHeadroomFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(outHzFinal(OutputGABARIndexs1D), GABAHeadroomFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Out Hz'); ylabel('GABA headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    print([directory,GABA_prep,'components_vs_headroom'],'-dpng');
    print([directory,GABA_prep,'components_vs_headroom'],'-depsc');
    
    %% GABA time to headroom
    figure(GABATimeToHeadroomFigNum);
    
    subplot(2,7,1);    
    scatter(AMPAWeights1D, GABATimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(AMPAWeights1D, GABATimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('AMPA'); ylabel('GABA time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    subplot(2,7,8);    
    scatter(GABAWeights1D, GABATimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABAWeights1D, GABATimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA'); ylabel('GABA time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,2);
    scatter(glutamateHzFinal(glutamatePreIndexs1D), GABATimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateHzFinal(glutamatePreIndexs1D), GABATimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Glutamate Hz'); ylabel('GABA time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,9);
    scatter(GABAHzFinal(GABAPreIndexs1D), GABATimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABAHzFinal(GABAPreIndexs1D), GABATimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA Hz'); ylabel('GABA time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,3);
    scatter(glutamateCB1RFinal, GABATimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateCB1RFinal, GABATimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Glutamate CB1R'); ylabel('GABA time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,10);
    scatter(GABACB1RFinal, GABATimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABACB1RFinal, GABATimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA CB1R'); ylabel('GABA time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,4);
    scatter(glutamateAvailableNeurotransmitterFinal, GABATimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateAvailableNeurotransmitterFinal, GABATimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Available glutamate'); ylabel('GABA time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,11);
    scatter(GABAAvailableNeurotransmitterFinal, GABATimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABAAvailableNeurotransmitterFinal, GABATimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Available GABA'); ylabel('GABA time to headroom');
    set([gca; findall(gca, 'Type','texglutamatePreIndexs1Dt')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,5);
    scatter(spineCaFinal, GABATimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(spineCaFinal, GABATimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Spine Ca2+'); ylabel('GABA time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,12);
    scatter(dendriteCaFinal, GABATimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(dendriteCaFinal, GABATimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Dendrite Ca2+'); ylabel('GABA time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,6);
    scatter(glutamateCleftAstrocyteeCBFinal, GABATimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateCleftAstrocyteeCBFinal, GABATimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Glutamate cleft eCB'); ylabel('GABA time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,13);
    scatter(GABACleftAstrocyteeCBFinal, GABATimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABACleftAstrocyteeCBFinal, GABATimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA cleft eCB'); ylabel('GABA time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
        
    subplot(2,7,[7,14]);
    scatter(outHzFinal(OutputGABARIndexs1D), GABATimeToHeadroom,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(outHzFinal(OutputGABARIndexs1D), GABATimeToHeadroom);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Out Hz'); ylabel('GABA time to headroom');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    print([directory,GABA_prep,'components_vs_timeToHeadroom'],'-dpng');
    print([directory,GABA_prep,'components_vs_timeToHeadroom'],'-depsc');
    
    %% Component scatters
    [~, figNum] = newFigure(figNum, true); 
    
    subplot(2,2,1);
    scatter(AMPAWeights1D, glutamateHzFinal(glutamatePreIndexs1D),'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(AMPAWeights1D, glutamateHzFinal(glutamatePreIndexs1D));
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('AMPA weight'); ylabel('Glutamate Hz');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    subplot(2,2,2);
    scatter(GABAWeights1D, glutamateHzFinal(glutamatePreIndexs1D),'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABAWeights1D, glutamateHzFinal(glutamatePreIndexs1D));
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA weight'); ylabel('Glutamate Hz');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    subplot(2,2,3);
    scatter(AMPAWeights1D, GABAHzFinal(GABAPreIndexs1D),'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(AMPAWeights1D, GABAHzFinal(GABAPreIndexs1D));
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('AMPA weight'); ylabel('GABA Hz');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    subplot(2,2,4);
    scatter(GABAWeights1D, GABAHzFinal(GABAPreIndexs1D),'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABAWeights1D, GABAHzFinal(GABAPreIndexs1D));
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA weight'); ylabel('GABA Hz');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    print([directory,'weights_vs_Hz'],'-dpng');
    print([directory,'weights_vs_Hz'],'-depsc');
    
    %% Component scatters
    [~, figNum] = newFigure(figNum, true); 
    
    subplot(2,2,1);
    scatter(glutamateCleftAstrocyteeCBFinal, glutamateCB1RFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateCleftAstrocyteeCBFinal, glutamateCB1RFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Glutamate eCB'); ylabel('Glutamate CB1R');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    subplot(2,2,2);
    scatter(glutamateCB1RFinal, AMPAWeights1D,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(glutamateCB1RFinal, AMPAWeights1D);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Glutamate CB1R'); ylabel('AMPA weight');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared; 
    
    subplot(2,2,3);
    scatter(GABACleftAstrocyteeCBFinal, GABACB1RFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABACleftAstrocyteeCBFinal, GABACB1RFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA eCB'); ylabel('GABA CB1R');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    subplot(2,2,4);
    scatter(GABACB1RFinal, GABAWeights1D,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(GABACB1RFinal, GABAWeights1D);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('GABA CB1R'); ylabel('GABA weight');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    print([directory,'eCB_vs_CB1R_vs_W'],'-dpng');
    print([directory,'eCB_vs_CB1R_vs_W'],'-depsc');
    
    %% Component scatters 2
    [~, figNum] = newFigure(figNum, true); 
    
    subplot(2,2,1);
    scatter(spineCaFinal, glutamateCleftAstrocyteeCBFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(spineCaFinal, glutamateCleftAstrocyteeCBFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Spine Ca2+'); ylabel('Glutamate eCB');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    subplot(2,2,2);
    scatter(spineCaFinal, GABACleftAstrocyteeCBFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(spineCaFinal, GABACleftAstrocyteeCBFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Spine Ca2+'); ylabel('GABA eCB');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared; 
    
    subplot(2,2,3);
    scatter(dendriteCaFinal, glutamateCleftAstrocyteeCBFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(dendriteCaFinal, glutamateCleftAstrocyteeCBFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Dendrite Ca2+'); ylabel('Glutamate eCB');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    subplot(2,2,4);
    scatter(dendriteCaFinal, GABACleftAstrocyteeCBFinal,'.');
    h=lsline;
    set(h,'color','r');
    R=corrcoef(dendriteCaFinal, GABACleftAstrocyteeCBFinal);
    R_squared=R(2)^2;
    title(['(R squared= ', num2str(R_squared), ')']);
    xlabel('Dendrite Ca2+'); ylabel('GABA eCB');
    set([gca; findall(gca, 'Type','text')], 'FontSize', fontsize);
    clear h R R_squared;
    
    print([directory,'Ca_vs_eCB'],'-dpng');
    print([directory,'Ca_vs_eCB'],'-depsc');
    
    %%
end
function [glutamateFluxPerturbation, glutamateFluxPerturbationMean, ...
    glutamateFluxPerturbationMeanBinned, glutamateFluxPerturbationP] = ...
    glutamateFluxVsPerturbation(cleftGlutamate, afterTimeToHeadroom, ...
    diffPerturbation, headroomConverged, p, HzMax, AMPAMax, perturbationType, ...
    titleStr)
    subplot(10,2,1:2:20); hold on;
    glutamateFluxPerturbation = zeros(1,size(cleftGlutamate,2));
    for i=1:size(cleftGlutamate,2)
        glutamateFluxPerturbation(i) = ...
            sum(cleftGlutamate(1:afterTimeToHeadroom(i),i));
    end
    scatter(diffPerturbation(headroomConverged), ...
        glutamateFluxPerturbation(headroomConverged),'.');
    if (p==1)
        Range = -HzMax+5:2:HzMax-5;
        offset = 5;
    elseif (p==2)
        Range = -AMPAMax+0.1:0.05:AMPAMax-0.1;
        offset = 0.1;
    end
    glutamateFluxPerturbationMean = zeros(1,numel(Range));
    tempDiffPerturbation = diffPerturbation;
    i = 1;
    for R=Range
        temp = tempDiffPerturbation >= R-offset & ...
            tempDiffPerturbation < R+offset;
        glutamateFluxPerturbationMean(i) = ...
            mean(glutamateFluxPerturbation(temp));
        i=i+1;
    end
    plot(Range, glutamateFluxPerturbationMean, 'LineWidth', 2);
    if (p==1)
        xlim([-HzMax HzMax]);
    elseif (p==2)
        xlim([-AMPAMax AMPAMax]);
    end
    title([perturbationType, ' perturbation vs glutamate flux', titleStr]);
    xlabel([perturbationType,' perturbation']);
    ylabel('glutamate flux');
    clear Range offset;

    subplot(10,2,2:2:17);
    if (p==1)
        Range = -140:20:140;
        offset = 10;
    elseif (p==2)
        Range = -1.4:0.2:1.4;
        offset = 0.1;
    end
    glutamateFluxPerturbationMeanBinned = zeros(1,numel(Range));
    tempDiffPerturbation = diffPerturbation;
    glutamateFluxPerturbationP = ones(1,numel(Range));
    % Find little change case for stat tests
    littleChange = glutamateFluxPerturbation(tempDiffPerturbation >= -offset & ...
        tempDiffPerturbation < offset);
    for i=1:numel(Range)
        % Find mean and stat to little change case
        temp = tempDiffPerturbation >= Range(i)-offset & ...
            tempDiffPerturbation < Range(i)+offset;
        glutamateFluxPerturbationMeanBinned(i) = ...
            mean(glutamateFluxPerturbation(temp));
        if (numel(glutamateFluxPerturbation(temp)) > 0 && ...
            numel(littleChange > 0))
            [glutamateFluxPerturbationP(i)] = ...
                ranksum(glutamateFluxPerturbation(temp), ...
                littleChange);
        end
    end
    bar(Range, glutamateFluxPerturbationMeanBinned);
    title([perturbationType,' perturbation vs glutamate flux', titleStr]);
    xlabel([perturbationType,' perturbation']);
    if (p==1)
        xlim([-HzMax+10 HzMax-10]);
        xticks([-140:40:-20,0,20:40:140]);
    elseif (p==2)
        xlim([-AMPAMax AMPAMax]);
        xticks(-1.4:0.2:1.4);
    end
    ylabel('glutamate flux');

    subplot(10,2,20);
    plot(Range, glutamateFluxPerturbationP<0.005,'*');
    set(gca,'XTick',[]); set(gca,'YTick',[]);
    if (p==1)
        xlim([-HzMax+10 HzMax-10]);
    elseif (p==2)
        xlim([-AMPAMax AMPAMax]);
    end
    ylim([0.5 1.5]);
end
