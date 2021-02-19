% xCoMP Configuration
% read configuration struct from system_cfg.xml
% tmpStruct = xml2struct(configureXML);
% systemCfg = pruneXMLstruct(tmpStruct.SystemCfg);

rng(systemCfg.SLS_Top_Cfg.RAND_SEED);
num_drops=systemCfg.SLS_Top_Cfg.NUM_DROPS;
num_iter=systemCfg.SLS_Top_Cfg.NUM_FADING_REALIZATIONS;

deployment_model=systemCfg.NetworkLayout.LayoutModel;

%if plot stats
isStats_plot=systemCfg.SLS_Top_Cfg.isStatsPlotting;
fig_stat=30;%*****FIXME: to be replaced later
color='b';%****b
isSave=systemCfg.SLS_Top_Cfg.isSave;%****
Quadriga_used = systemCfg.ChannelSettings.Quadriga_used;


%Drop conifigurations
BLO=systemCfg.NetworkLayout.BLO;
LO=systemCfg.NetworkLayout.DroppingModel;%"RegRndm_X";%"RegRndm";%;%"RndmRndm";%"HexRndm";
% FIXME: need to modify the logics below
if BLO
    layout_len=systemCfg.NetworkLayout.Layout_Length;
    layout_wid=systemCfg.NetworkLayout.Layout_Width;
    Lo=systemCfg.NetworkLayout.Lo;
    Wo=systemCfg.NetworkLayout.Wo;
    Nx=systemCfg.NetworkLayout.RelayNodes.Nx; %better be even
else
    layout_len=systemCfg.NetworkLayout.Layout_Length;
    layout_wid=systemCfg.NetworkLayout.Layout_Width;
    Lo=systemCfg.NetworkLayout.Lo;
    Wo=systemCfg.NetworkLayout.Wo;
    Nx=systemCfg.NetworkLayout.RelayNodes.Nx; %better be even
end
height_struct=struct('gnb',systemCfg.NetworkLayout.Height_gnb,'ue',systemCfg.NetworkLayout.Height_ue,'sc',systemCfg.NetworkLayout.RelayNodes.Height_x); %in Cmeters above UE height
num_gNB_rows=systemCfg.gNB.NumGNb_Rows;
offset_x=0;%12.5; FIXME: to load from xml
isLOPlot=systemCfg.SLS_Top_Cfg.isLOPlotting;
fig_idx=20;


Nofl=systemCfg.NetworkLayout.RelayNodes.Nofl;
Nax=systemCfg.NetworkLayout.RelayNodes.Nax;


LO_dim=struct('LO_len',layout_len,'LO_wid', layout_wid,'Ht',height_struct,'Nrows',num_gNB_rows,'Nx',Nx,'N_offload',Nofl,'Lo',Lo,'Wo',Wo);
rrue_type=systemCfg.NetworkLayout.RelayNodes.RRUE_Type;%"dim_shp_rept";%'dim_shp';%'snr_shp';

%%% scenario related
scenario = systemCfg.NetworkLayout.Scenario;
bler_target = systemCfg.gNB.BLER_Target;

% Channel aspects
Nc_store=systemCfg.ChannelSettings.Nc_store;
%below is to store channel realizations and run new sims on the same H
storeChan=struct('isReadChannel',0,'isSaveChannel',0,'maxCol',Nc_store);

if and(storeChan.isReadChannel,storeChan.isSaveChannel)
    assert;
end

fc=systemCfg.ChannelSettings.CarrierFrequency/1e9; %in GHz
fadingModel=systemCfg.ChannelSettings.FadingModel;
%%% Read channel setting related parameters, which depends on if BW is indicated from RB perspective or Hz perspective
if systemCfg.ChannelSettings.RBEnable == 1
    BW_RB=systemCfg.ChannelSettings.SystemBandwidth;  %%% in unit of RBs
    fade_gran_RB=systemCfg.ChannelSettings.FadingGranularity; %%% in unit of RBs
    BW = 12 * BW_RB * systemCfg.ChannelSettings.SubCarrierSpacing / 1e6; %in MHz
    fade_gran = 12 * fade_gran_RB * systemCfg.ChannelSettings.SubCarrierSpacing / 1e6; %in MHz
else
    BW=systemCfg.ChannelSettings.SystemBandwidth/1e6;%in MHz
    fade_gran=systemCfg.ChannelSettings.FadingGranularity/1e6; %BW of a fade in MHz
end
del_f=systemCfg.ChannelSettings.SubCarrierSpacing/1e3; %inter carrier spacing in KHz
Nfreq=ceil(BW/fade_gran);   %%% Wanlu: replace floor by ceil
Krician=systemCfg.ChannelSettings.K_Rician;%10^(10/10);

%antenna configurations
%tput factor
Sys_Tput_factor = 1200*12*0.9/1000 * (BW/20); %%% Wanlu: what does BW/20 mean? need to be revised now?

if BLO
    %%%%%%%BLO%%%%%
    Nu = systemCfg.NetworkLayout.Num_UE;
    Nr=systemCfg.UE.NumUE_Antenna; %per UE
    Nr_ob=systemCfg.UE.NumUE_Antenna_Observed;
    L=systemCfg.UE.NumMIMOLayers;
    Nt=systemCfg.gNB.NumGNb_Antenna; %per AP
    Nb = systemCfg.NetworkLayout.Num_gNB;
    if (systemCfg.NetworkLayout.isGNB_SquareLayout)
        num_gNB_rows=floor(sqrt(Nb));
    end
else
    Nu = systemCfg.NetworkLayout.Num_UE;
    Nr=systemCfg.UE.NumUE_Antenna; %per UE
    Nr_ob=systemCfg.UE.NumUE_Antenna_Observed;
    L=systemCfg.UE.NumMIMOLayers;
    Nt=systemCfg.gNB.NumGNb_Antenna; %per AP
    Nb = systemCfg.NetworkLayout.Num_gNB;
    if (systemCfg.NetworkLayout.isGNB_SquareLayout)
        num_gNB_rows=floor(sqrt(Nb));
    end
end

precoding_type=systemCfg.gNB.BBU.BeamForming.BF_Algo;%"BDiag";%"trunc_SLR";%"orth-orth";%"RZF";%"orth-orth";%"ZF";%"orth-reg";"random";******
lambda=systemCfg.gNB.BBU.BeamForming.RZF_lambda;%only used for RZF, and orth-reg ******
txMode=systemCfg.gNB.BBU.BeamForming.txMode;%"mMIMO";%"CoMP";%"mMIMO";%%%;%C; ****

%link considerations
if BLO
    Ptx=systemCfg.gNB.PTx;%in dBm  ******
    Ptxx=systemCfg.gNB.PTxRelay;
else
    Ptx=systemCfg.gNB.PTx;%in dBm  ******
    Ptxx=systemCfg.gNB.PTxRelay;
end
PtxUE=systemCfg.UE.PTxUE;

%FIXME: to update Tx/Rx gain and noise figure, etc.
if isfield(systemCfg.gNB.AntennaPattern, 'GTx_gNB_dB')
Gt=systemCfg.gNB.AntennaPattern.GTx_gNB_dB;% in dB
else
    Gt=systemCfg.gNB.GTx_gNB_dB;% in dB
end
Gr=systemCfg.UE.GRx_UE_dB;%in dB

if isfield(systemCfg.gNB.AntennaPattern, 'GTx_Relay_dB')
Gtx=systemCfg.gNB.AntennaPattern.GTx_Relay_dB;% in dB
else
    Gtx=systemCfg.gNB.GTx_Relay_dB;% in dB
end


noise_figure=systemCfg.UE.NoiseFigure; %in dB
MCL_struct=struct('gnb',systemCfg.NetworkLayout.MCL,'sc',systemCfg.NetworkLayout.RelayNodes.MCLx);%in dBm, min coupling loss
MCL=systemCfg.NetworkLayout.MCL;

% %%% old codes
% PtxG=Ptx+Gt+Gr;%-30;%30 to convert to dBp from dBm
% PtxGx=Ptxx+Gtx+Gr;
% PtxUEG=PtxUE+Gt+Gr;
% %%%%%%%%%%%%%%

%%% new codes
PtxG=Ptx+Gr;%-30;%30 to convert to dBp from dBm
PtxGx=Ptxx+Gr;
PtxUEG=PtxUE+Gr;
%%%%%%%%%%%%%%

bypass_shadowing=systemCfg.NetworkLayout.isBypassShadowing;
noise_floor=-174+10*log10(BW*1e6)+noise_figure; %noise floor in dBm
SNR_cap=systemCfg.UE.SNRdB_Cap;%cap on attained SNR
gap2cap=systemCfg.UE.gap2Cap;%in dB
gap2cap_lin=10^(gap2cap/10);
Rx_coloring=1;

if (strcmp(systemCfg.UE.Demapper, 'MMSE_IRC'))
    isPerantennaSNR=1;
elseif (strcmp(systemCfg.UE.Demapper, 'MMSE_SIC'))
    isPerantennaSNR=0;
else
    error('Invalid Demapper Mode\n');
end

%channel error model
if (strcmp(systemCfg.gNB.RRU.Calib_Err, 'ON'))
    chanPerturb=1;%******
else
    chanPerturb=0;%******
end
isChan_err=systemCfg.gNB.RRU.isChan_Err;%error happens to all channel entries not only for diagonal for calib
phase_err=systemCfg.gNB.RRU.Phase_Err; %error in degree will UD[-phase_err,phase_err]--> mean=0, variance=(phase_err^2)/3, and rms=phase_err/sqrt(3)
mag_err=10^(systemCfg.gNB.RRU.Mag_Err_dB/10);%the variance of one components of the compelx Gaussian error, i.e. error=sqrt(mag_err)(a+jb), where a & b are N(0,1)
%Magnitude is Rayleigh(mag_err), that is mean=mag_errxsqrt(pi/2) and
%var=((4-pi)/2)mag_err^2
flag_rf_phase_error = 0;

%AS related
overdim_factor=systemCfg.gNB.BBU.Cell_Association.OverDimensionFactor;
isUE_Assoc_Enabled = systemCfg.gNB.BBU.Cell_Association.isUE_Assoc_Enabled;
AS_metric=systemCfg.gNB.BBU.Cell_Association.AS_Metric;%'ASdimMetric';%"ASCapMetric";%"SLR";%"ASCapMetric";%"SLR";%"ASCapMetric";%"ASdimMetric";%'ASlinqMetric';%"ASCapMetric";%"ASdimMetric";%"ASCapMetric";%'ASlinqMetric';%;%'ASlinqMetric',"ASCapMetric"
if BLO
    AS_size_max=systemCfg.gNB.BBU.Cell_Association.AS_Size_Max;
else
    AS_size_max=systemCfg.gNB.BBU.Cell_Association.AS_Size_Max;
end
O_complexity=systemCfg.gNB.BBU.Cell_Association.O_complexity;
lkg_exponent=systemCfg.gNB.BBU.Cell_Association.lkg_exponent;
max_lkg_allowed=systemCfg.gNB.BBU.Cell_Association.max_lkg_allowed;% in dBc tp slr_max
if ~isfield(systemCfg.gNB.BBU.Cell_Association, 'clustering_type') && ...
        isfield(systemCfg.gNB.BBU.Cell_Association, 'flag_single_cluster')
    if systemCfg.gNB.BBU.Cell_Association.flag_single_cluster == 1
        systemCfg.gNB.BBU.Cell_Association.clustering_type = 'single';
    else
        assert(false);
    end
end
clustering_type = systemCfg.gNB.BBU.Cell_Association.clustering_type;

dims=struct('Nt',Nt,'Nr',Nr,'Nr_ob',Nr_ob,'Nu',Nu,'Nb',Nb,'Nf',ceil(BW/fade_gran),'Ncm',storeChan.maxCol);   %%% Wanlu: replace floor by ceil
Lo=LO_dim.Lo;
Wo=LO_dim.Wo;
Nx=LO_dim.Nx;

% some sanity checks
if (Nt*Nb<Nr_ob*Nu) && ~systemCfg.Scheduler.Enabled
    assert;
end


%UL
UL_Clustering=systemCfg.UL.UL_Clustering;
isUL=  systemCfg.SLS_Top_Cfg.isUL;
filter_type_UL=systemCfg.UE.ULfilterType;
%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%Hetro Deployment (Macro+hotspots)

Lx=systemCfg.HetNet.hotspot_len;
Wx=systemCfg.HetNet.hotspot_wid;
numx_rows=systemCfg.HetNet.hotspot_numRows;
isTrackMacroInterf=systemCfg.HetNet.hotspot_isTrackMacroInterf;
fig_stat_hetro=systemCfg.HetNet.hotspot_figStat_num;
isfixed_sc_macro_LO=systemCfg.HetNet.hotspot_isFixedScMacroLO;
ueDensity_sc=systemCfg.HetNet.hotspot_ueDensity;
lenC=systemCfg.HetNet.hotspot_numHotspots;
hetro_reuse1=systemCfg.HetNet.hotspot_reuse1;
boundary_mx=100+Lx;
boundary_my=100+Wx;


% All parameters related to mmWave simulations
if ~isfield(systemCfg.ChannelSettings, 'Quadriga_used')
    systemCfg.ChannelSettings.Quadriga_used =0;
end

if systemCfg.ChannelSettings.Quadriga_used
    Nt2 = systemCfg.gNB.NumGNb_Antenna_AZ;  % added by Ali to account for the number of ant. in other dimsion of the panel. If linear array, then Nt2 =0;
    minDistance = systemCfg.NetworkLayout.VOIDLayout.minDistanceBetweenGroups;
    groupRaduis = systemCfg.NetworkLayout.VOIDLayout.groupRaduis;
    num_group   =  systemCfg.NetworkLayout.VOIDLayout.num_group;
    UE_Pattern   =  systemCfg.NetworkLayout.VOIDLayout.UE_Pattern;
    deply_SC = systemCfg.NetworkLayout.VOIDLayout.deployment;
    Quadriga_used = systemCfg.ChannelSettings.Quadriga_used;
    Blkg_dur = systemCfg.NetworkLayout.VOIDLayout.Blkg_dur;
    gradualBlkg = systemCfg.NetworkLayout.VOIDLayout.gradualBlkg;
    %void
    SNR_measuremnt = systemCfg.NetworkLayout.VOIDLayout.SNR_measuremnt;
    SLRBeamCoordEnabled = systemCfg.NetworkLayout.VOIDLayout.SLRBeamCoordEnabled;
    QuarterBlokageEnabled = systemCfg.NetworkLayout.VOIDLayout.QuarterBlokageEnabled;
    beamAging    =   systemCfg.NetworkLayout.VOIDLayout.beamAging;
    TopBeams    =   systemCfg.NetworkLayout.VOIDLayout.TopBeams;
    Slicing    =   systemCfg.NetworkLayout.VOIDLayout.Slicing;
    reTXdiffAP    =   systemCfg.NetworkLayout.VOIDLayout.retransmission.reTXdifferentAP;
    TopBeamsDiffQuadrants  =  systemCfg.NetworkLayout.VOIDLayout.retransmission.TopBeamsDiffQuadrants;
    reTXsuccessiveFailure    =   systemCfg.NetworkLayout.VOIDLayout.retransmission.reTXsuccessiveFailure;
    BAfactor    =   systemCfg.NetworkLayout.VOIDLayout.BAfactor;
    UE_orientation = systemCfg.NetworkLayout.VOIDLayout.UE_orientation;
    TrackCoordinatesFromFile =  systemCfg.NetworkLayout.VOIDLayout.TrackCoordinatesFromFile;
    TrackSectionNumber  =  systemCfg.NetworkLayout.VOIDLayout.TrackSectionNumber;
    secNumber         = systemCfg.NetworkLayout.VOIDLayout.TrackSectionNumber;
    st80211ad = systemCfg.techStandard.st80211ad;
    TargetError = systemCfg.techStandard.TargetError;
    st80211NR = systemCfg.techStandard.st80211NR;
    TimeSlot = systemCfg.techStandard.TimeSlot;
    Fixed_codebook = systemCfg.gNB.BBU.BeamForming.BF_FixedCodeBook;
    FixedCodeBook_type=systemCfg.gNB.BBU.BeamForming.BF_FixedCodeBook;
    if systemCfg.UE.UEHeight == 'Fixed'  %
        coordu_z = 1.65*ones(1, Nu);
    elseif systemCfg.UE.UEHeight == 'Variable'
        coordu_z = [1.5*ones(1, 0.75*Nu) 1*ones(1, 0.25*Nu)];
    end
end


%%% Parameters related to Hexagon cell deployment
if systemCfg.NetworkLayout.Num_gNB==7
    hotspot_numHotspots = [systemCfg.HetNet.hotspot_numHotspots1, systemCfg.HetNet.hotspot_numHotspots2, systemCfg.HetNet.hotspot_numHotspots3, ...
        systemCfg.HetNet.hotspot_numHotspots4, systemCfg.HetNet.hotspot_numHotspots5, systemCfg.HetNet.hotspot_numHotspots6, systemCfg.HetNet.hotspot_numHotspots7];
    
    hotspot_distmacro = [systemCfg.HetNet.hotspot_distmacro1, systemCfg.HetNet.hotspot_distmacro2, systemCfg.HetNet.hotspot_distmacro3, ...
        systemCfg.HetNet.hotspot_distmacro4, systemCfg.HetNet.hotspot_distmacro5, systemCfg.HetNet.hotspot_distmacro6, systemCfg.HetNet.hotspot_distmacro7];
elseif systemCfg.NetworkLayout.Num_gNB==19
    hotspot_numHotspots = [systemCfg.HetNet.hotspot_numHotspots1, systemCfg.HetNet.hotspot_numHotspots2, systemCfg.HetNet.hotspot_numHotspots3, ...
        systemCfg.HetNet.hotspot_numHotspots4, systemCfg.HetNet.hotspot_numHotspots5, systemCfg.HetNet.hotspot_numHotspots6, ...
        systemCfg.HetNet.hotspot_numHotspots7, systemCfg.HetNet.hotspot_numHotspots8, systemCfg.HetNet.hotspot_numHotspots9, ...
        systemCfg.HetNet.hotspot_numHotspots10, systemCfg.HetNet.hotspot_numHotspots11, systemCfg.HetNet.hotspot_numHotspots12, ...
        systemCfg.HetNet.hotspot_numHotspots13, systemCfg.HetNet.hotspot_numHotspots14, systemCfg.HetNet.hotspot_numHotspots15, ...
        systemCfg.HetNet.hotspot_numHotspots16, systemCfg.HetNet.hotspot_numHotspots17, systemCfg.HetNet.hotspot_numHotspots18, systemCfg.HetNet.hotspot_numHotspots19];
    
    hotspot_distmacro = [systemCfg.HetNet.hotspot_distmacro1, systemCfg.HetNet.hotspot_distmacro2, systemCfg.HetNet.hotspot_distmacro3, ...
        systemCfg.HetNet.hotspot_distmacro4, systemCfg.HetNet.hotspot_distmacro5, systemCfg.HetNet.hotspot_distmacro6, ...
        systemCfg.HetNet.hotspot_distmacro7, systemCfg.HetNet.hotspot_distmacro8, systemCfg.HetNet.hotspot_distmacro9, ...
        systemCfg.HetNet.hotspot_distmacro10, systemCfg.HetNet.hotspot_distmacro11, systemCfg.HetNet.hotspot_distmacro12, ...
        systemCfg.HetNet.hotspot_distmacro13, systemCfg.HetNet.hotspot_distmacro14, systemCfg.HetNet.hotspot_distmacro15, ...
        systemCfg.HetNet.hotspot_distmacro16, systemCfg.HetNet.hotspot_distmacro17, systemCfg.HetNet.hotspot_distmacro18, systemCfg.HetNet.hotspot_distmacro19];
elseif systemCfg.NetworkLayout.Num_gNB==1
    hotspot_numHotspots = max(systemCfg.HetNet.hotspot_numHotspots1, systemCfg.HetNet.hotspot_numHotspots);
    hotspot_distmacro = systemCfg.HetNet.hotspot_distmacro1;
else
    hotspot_numHotspots = 0;
    hotspot_distmacro = 0;
%     assert(1<0, 'Wrong cell deployment.');
end
hotspot_dist_limit = systemCfg.HetNet.hotspot_margin +sqrt((Lx)^2+(Wx)^2) ;


SC_outdoor=struct('Lx',Lx,'Wx',Wx,'numx_rows',numx_rows,'boundary_mx',boundary_mx,'boundary_my',boundary_my,'reuse1',hetro_reuse1,...
    'isTrackMacroInterf',isTrackMacroInterf,'isfixed_sc_macro_LO',isfixed_sc_macro_LO,'ueDensity_sc',ueDensity_sc, 'hotspot_distmacro', hotspot_distmacro, 'hotspot_numHotspots', hotspot_numHotspots, ...
    'hotspot_place', systemCfg.HetNet.hotspot_place, 'hotspot_dist_limit', hotspot_dist_limit, 'lenC', lenC, 'Wx_ue', systemCfg.HetNet.hotspot_ueDens_wid, 'Lx_ue', systemCfg.HetNet.hotspot_ueDens_len, 'Num_sector', systemCfg.HetNet.Num_sector,...
    'Triangle_shape', systemCfg.HetNet.Triangle_shape);


systemCfg.NetworkLayout.SC_outdoor=SC_outdoor;
systemCfg.HetNet.hotspot_numHotspots = hotspot_numHotspots;