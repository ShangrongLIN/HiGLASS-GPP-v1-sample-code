clear all;
% the origin path should be here please change with your demand
originPATH='H:\HIGLASS_test\HIGLASS_GPP_MATLAB_TEST';

outputPATH=strcat(originPATH,'\OUTPUTdata');
mkdir(outputPATH);
TEMdata=imread(strcat(originPATH,'\CLIMATEdata\ERA_TEMP_month_degreeC_ori_16to20.tif'));
VPDdata=imread(strcat(originPATH,'\CLIMATEdata\ERA_VPD_month_hPa_ori_16to20.tif'));
PAdata=imread(strcat(originPATH,'\CLIMATEdata\MERRA_meanPA_month_Pa_ori_16to20.tif'));
PARdata=imread(strcat(originPATH,'\CLIMATEdata\MERRA_sumPAR_month_MJ_ori_16to20.tif'));





datapath_NDVI=strcat(originPATH,'\NDVI');
datapath_LUCC=strcat(originPATH,'\LUCC2');
datapath_LUCC_geo_trans=strcat(originPATH,'\LUCC_geo');
datapath_MAIZE=strcat(originPATH,'\Maize');

folder_NDVI=dir(datapath_NDVI);
folder_LUCC=dir(datapath_LUCC);

%runTAG=11;  %testTAGhere, represent for the path less than 1
for foldernum=3:length(folder_NDVI)
       subfolderpath_NDVI= strcat(datapath_NDVI,'\',folder_NDVI(foldernum).name);
        
       subfolder_NDVI=dir(subfolderpath_NDVI);
       %subfolderpath_LUCC=strcat(datapath_LUCC,'\',folder_LUCC(foldernum).name);
       %subfolder_LUCC=dir(subfolderpath_LUCC);
       %F:\HIGLASS_PRODUCT\run_data\run_region_test_WHU\NDVI_month_134039134040_201601202009\hh26vv06\h23v08
       
       %% merge the NDVI data here, since a lot of data missed
       for subfoldernum=3:3
           %  read the data and composite them into a mat with 927*927*60
           %  month
              %sub_subfolderpath= strcat(subfolderpath_NDVI,'\',subfolder_NDVI(subfoldernum).name);
              % a lot of data are missed here, so we just keep the data with  available data
              [saveTAG]=findlistTAG(subfolder_NDVI);
              NDVImat=int16(zeros(927,927,60));
              QCmat=uint8(zeros(927,927,60));
              for tag=1:length(saveTAG)
                  if saveTAG(tag,1)>0
                     filename_general=subfolder_NDVI(saveTAG(tag,1)).name;
                     filename_subset=filename_general(6:19);
                     NDVI_filename=strcat(subfolderpath_NDVI,'\',filename_general,'\',filename_general,'_ndvi30.tif');
                     QC_filename=strcat(subfolderpath_NDVI,'\',filename_general,'\',filename_general,'_ndvi30_QA.tif');
                     NDVImat(1:927,1:927,tag)=imread(NDVI_filename);
                     QCmat(1:927,1:927,tag)=imread(QC_filename);
                  else  
                     NDVImat(1:927,1:927,tag)=int16(ones(927,927)*-9999);   % define a standard size
                     QCmat(1:927,1:927,tag)=uint8(ones(927,927)*2);   % define a standard size
                  end
              end
             % so start to rebuilt the data here
       end

        
        %% rebuilt the NDVI data with QC
        [rebuiltNDVI,rebuiltQC]=rebuiltNDVI_fun(NDVImat,QCmat);
        
        clear NDVImat QCmat       

       %% find the LUCC data
        LUCC_datapath_temp=dir(strcat(datapath_LUCC,'\GLC_',filename_subset,'_FCS30_2020_*'));
        LUCC_datapath=strcat(datapath_LUCC,'\',LUCC_datapath_temp(1).name);
        [LUCCdata,geoR]=geotiffread(LUCC_datapath);
        MaizeFilename=strcat(datapath_MAIZE,'\maize-',filename_subset,'_FCS30-2020.tif');
        [Maizedata,~]=geotiffread(MaizeFilename);
        [LUCC_EC_LUE] = changeLUCC(Maizedata,LUCCdata);
        %% get the climate data
        % first we need to get the longitude and latitude of each pixel to get the size
        [~,geoR_latlon]=geotiffread(strcat(datapath_LUCC_geo_trans,'\',LUCC_datapath_temp(1).name,'_GEO.tif'));
        
         image_region_limit_lon=geoR_latlon.LongitudeLimits;
         image_region_limit_lat=geoR_latlon.LatitudeLimits;
         

         [T_timeseries] = get_T_region_openDATA(TEMdata,image_region_limit_lat,image_region_limit_lon);
         [VPD_timeseries] = get_VPD_region_openDATA(VPDdata,image_region_limit_lat,image_region_limit_lon);
         [PA_timeseries] = get_PA_region_openDATA(PAdata,image_region_limit_lat,image_region_limit_lon);
         [PAR_timeseries] = get_PAR_region_openDATA(PARdata,image_region_limit_lat,image_region_limit_lon);
        %% run EC-LUE
        
        
         
        GPP_time_series_mat=higlass_EC_LUE_function(LUCC_EC_LUE,rebuiltNDVI,T_timeseries,VPD_timeseries,PAR_timeseries,PA_timeseries);
        outputGPPname=strcat(outputPATH,'\',filename_general(1:19),'_GPP_2016_2020_month.tif');
        outputQCname=strcat(outputPATH,'\',filename_general(1:19),'_GPP_QC_2016_2020_month.tif');
        info = geotiffinfo(LUCC_datapath);
        geotiffwrite(outputGPPname,GPP_time_series_mat,geoR,'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
        geotiffwrite(outputQCname,rebuiltQC,geoR,'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag);
        
        clear GPP_time_series_mat LUCC_EC_LUE LUCCdata  QCmat  rebuiltNDVI rebuiltQC LUCC_datapath
end



%% this is the main function for calculate EC-LUE based GPP


function [outputGPPmat]=higlass_EC_LUE_function(inputLUCC,NDVImonth_series,Tmonth_series,VPDmonth_series,PARmonth_series,PAmonth_series)
          CO2month_series=[402.46;403.41;403.55;404.78;404.42;404.59;404.23;404.40;404.85;405.22;405.73;405.33;406.05;405.82;406.06;406.38;406.38;406.69;407.00;407.29;407.16;407.21;407.34;407.71;408.85;408.68;408.05;407.62;408;408.60;408.57;409.19;409.32;409.55;410.20;409.95;410.74;411.12;410.65;410.71;411.42;411.77;411.64;412.21;412.37;412.11;412.47;412.66;413.31;413.50;413.20;413.61;413.87;414.21;414.29;414.81;415.13;414.86;415.11;414.94];
          input_spatial_size=927;
          outputGPPmat=uint16(zeros(input_spatial_size,input_spatial_size,60));
       for month=1:60
          Tmonth=resampleCLIMATEdata(Tmonth_series(month,1:9));
          VPDmonth=resampleCLIMATEdata(VPDmonth_series(month,1:9));
          PARmonth=resampleCLIMATEdata(PARmonth_series(month,1:9));
          PAmonth=resampleCLIMATEdata(PAmonth_series(month,1:9));          
          
          
          CO2month=CO2month_series(month,1);
          inputNDVI=NDVImonth_series(1:927,1:927,month);
          GPPpixel=uint16(zeros(input_spatial_size,input_spatial_size));
             
         % [inputNDVI,R2]=geotiffread(strcat(NDVIpath,'\',NDVIfile(month).name));
      if month==5 ||month==6 ||month==7 ||month==8 ||month==9 ||month==10 ||month==17 ||month==18 ||month==19||month==20||month==21 ||month==22 ||month==29 ||month==30 ||month==31||month==32||month==33 ||month==34||month==41 ||month==42 ||month==43||month==44||month==45 ||month==46||month==53 ||month==54 ||month==55||month==56||month==57 ||month==58
      % this is for maize algorithm
           for LUCCtype=1:10
               %%
               outputLUCCtemp=uint16(inputLUCC(1:927,1:927));
               switch LUCCtype
                   case 1  %maize
                    LUEmax=4.18;
                    theta=49.87;
                    VPD0=1.29;
                    Topt=25;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);
                    
                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
                   case 2  %C3 cro
                    LUEmax=2.93;
                    theta=49.87;
                    VPD0=1.29;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);
                    
                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;

                    %%
                   case 3 %EBF
                    LUEmax=2.78;
                    theta=29.81;
                    VPD0=0.5;
                    Topt=25;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);

                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
                    %%
                   case 4  %DBF
                    LUEmax=1.90;
                    theta=29.81;
                    VPD0=1.49;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);

                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
                    %%
                   case 5   %ENF
                    LUEmax=1.85;
                    theta=29.81;
                    VPD0=0.8;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);

                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
            %%        
                   case 6   %DNF
                    LUEmax=2.16;
                    theta=29.81;
                    VPD0=0.8;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);

                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
             %%       
                   case 7   %MF
                    LUEmax=2.33;
                    theta=51.95;
                    VPD0=1;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);

                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
           %%     
                   case 8   % SHR
                    LUEmax=2.05;
                    theta=36.72;
                    VPD0=1.39;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);

                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
            %%
                   case 9  % GRA
                    LUEmax=3.18;
                    theta=55.72;
                    VPD0=1.2;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);

                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
          %%          
                   case 10   %WET
                    LUEmax=3.18;
                    theta=55.72;
                    VPD0=1.20;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);
                    outputLUCCtemp=uint16(inputLUCC(1:927,1:927));
                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;             
               end
           end
      else
            for LUCCtype=2:10
               %%
               outputLUCCtemp=uint16(inputLUCC(1:927,1:927));
               switch LUCCtype
                    case 2  %C3 cro
                    LUEmax=2.93;
                    theta=49.87;
                    VPD0=1.29;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);
                    
                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;

                    %%
                   case 3 %EBF
                    LUEmax=2.78;
                    theta=29.81;
                    VPD0=0.5;
                    Topt=25;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);

                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
                    %%
                   case 4  %DBF
                    LUEmax=1.90;
                    theta=29.81;
                    VPD0=1.49;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);

                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
                    %%
                   case 5   %ENF
                    LUEmax=1.85;
                    theta=29.81;
                    VPD0=0.8;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);

                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
            %%        
                   case 6   %DNF
                    LUEmax=2.16;
                    theta=29.81;
                    VPD0=0.8;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);

                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
             %%       
                   case 7   %MF
                    LUEmax=2.33;
                    theta=51.95;
                    VPD0=1;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);

                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
           %%     
                   case 8   % SHR
                    LUEmax=2.05;
                    theta=36.72;
                    VPD0=1.39;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);

                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
            %%
                   case 9  % GRA
                    LUEmax=3.18;
                    theta=55.72;
                    VPD0=1.2;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);

                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;
          %%          
                   case 10   %WET
                    LUEmax=3.18;
                    theta=55.72;
                    VPD0=1.20;
                    Topt=20.33;
                    GPPtype_month=run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,Topt,theta,VPD0);
                    outputLUCCtemp=uint16(inputLUCC(1:927,1:927));
                    outputLUCCtemp(outputLUCCtemp~=LUCCtype)=0;
                    outputLUCCtemp(outputLUCCtemp==LUCCtype)=1;
                    GPPtype_month=outputLUCCtemp.*GPPtype_month;
                    GPPpixel=GPPpixel+GPPtype_month;                         
               end
           end  
           end
       
           outputGPPmat(:,:,month)=GPPpixel;   %GPP
           
       end
end

function [spatial_climate_data]=resampleCLIMATEdata(input_monthCLIMATE_DATA)
         x=1:3;
         y=1:3;
         temps=[input_monthCLIMATE_DATA(1,1)    input_monthCLIMATE_DATA(1,2)    input_monthCLIMATE_DATA(1,3)
                input_monthCLIMATE_DATA(1,4)    input_monthCLIMATE_DATA(1,5)    input_monthCLIMATE_DATA(1,6)
                input_monthCLIMATE_DATA(1,7)    input_monthCLIMATE_DATA(1,8)    input_monthCLIMATE_DATA(1,9)];
        xi=1:2/1853:3;
        yi=1:2/1853:3;    
        [xi,yi]=meshgrid(xi,yi);
        zi=interp2(x,y,temps,xi,yi,'linear');
        
        spatial_climate_data(1:927,1:927)=zi(465:1391,465:1391);
end


function [spatial_climate_data]=resampleCLIMATEdata_selfdefine(input_monthCLIMATE_DATA)

        xi=1:2/1853:3;
        yi=1:2/1853:3;    
        kx=((input_monthCLIMATE_DATA(1,3)+input_monthCLIMATE_DATA(1,9))/2-(input_monthCLIMATE_DATA(1,1)+input_monthCLIMATE_DATA(1,7))/2)/1854;
        ky=((input_monthCLIMATE_DATA(1,7)+input_monthCLIMATE_DATA(1,9))/2-(input_monthCLIMATE_DATA(1,1)+input_monthCLIMATE_DATA(1,3))/2)/1854;
        [xi,yi]=meshgrid(xi,yi);
        zi=input_monthCLIMATE_DATA(1,1)+xi*kx+yi*ky;

        spatial_climate_data(1:927,1:927)=zi(465:1391,465:1391);
end


function [GPPmonth] = run_EC_LUE_GPP(inputNDVI,Tmonth,VPDmonth,PARmonth,PAmonth,CO2month,LUEmax,topt,theta,VPD0)
PAR=PARmonth;
T=Tmonth;
VPD=VPDmonth;
patm=PAmonth;
CO2=CO2month;  % set a constant data here later need to change

kc25 = 39.97;      % Pa, assuming 25 deg C & 98.716 kPa
ko25 = 2.748e4;    % Pa, assuming 25 deg C & 98.716 kPa
dhac = 79430;      % J/mol
dhao = 36380;      % J/mol
kR   = 8.3145;     % J/mol/K
kco  = 2.09476e5;  % ppm, US Standard Atmosphere
Kc = kc25*exp(dhac*(T - 25.0)./(298.15*kR*(T + 273.15)));%ÉãÊÏ¶È
Ko = ko25*exp(dhao*(T - 25.0)./(298.15*kR*(T + 273.15)));
K  = Kc.*(1 + kco*(1.0e-6)*patm./Ko);
f1 = exp(-0.0227.*(T-25)); % the viscosity of water relative to its value at 25¡æ
epsilon=sqrt(356.51.*K./(1.6.*f1));
Ci=CO2.*epsilon./(epsilon+sqrt(VPD*100));%VPDÊäÈëpa  input data should be hPa
Cs=(Ci-theta)./(Ci+2.0.*theta);
%¼ÆËãts
tmin = 0.0;
tmax = 40.0;
%topt = 20.33;
ts = ((T-tmin).*(T-tmax))./((T-tmin).*(T-tmax)-(T-topt).*(T-topt));
ts = min(1, max(0, ts));
ws=VPD0./(VPD+VPD0);
ws = min(1, max(0, ws));
%     FPAR=1.24*NDVI_footprint-0.168;
FPAR=1.24*double(inputNDVI)/10000-0.168;
FPAR = min(1, max(0, FPAR));
GPPmonth_temp = LUEmax.*PAR.*FPAR.*Cs.* min(ts,ws);

GPPmonth=uint16(GPPmonth_temp*100);
end


 %% this function is for changing the LUCC
 function [outputLUCC] = changeLUCC(Maizedata,LUCC_map)
         %change the defined LUCC here
         %outputTEMP=inputLUCC
         LUCC_map=double(LUCC_map);
         LUCC_map(LUCC_map~=10 & LUCC_map~=11 & LUCC_map~=12  & LUCC_map~=20 & LUCC_map~=51 & LUCC_map~=52 & LUCC_map~=61 & LUCC_map~=62 & LUCC_map~=71 & LUCC_map~=72 & LUCC_map~=81 & LUCC_map~=82& LUCC_map~=91 & LUCC_map~=92 & LUCC_map~=120 & LUCC_map~=121 & LUCC_map~=122 & LUCC_map~=130 & LUCC_map~=180)=0;
         maize_map=double(Maizedata);
         CROmap_excludeMAIZE=(1-double(Maizedata));   % ä¸æ˜¯çŽ‰ç±³åˆ™å…¨éƒ¨æ˜¯1ï¼Œæ˜¯çŽ‰ç±³åˆ™å…¨éƒ¨æ˜¯0
         LUCC_map=LUCC_map.*CROmap_excludeMAIZE;       % mask æŽ‰çŽ‰ç±? 
         LUCC_map(LUCC_map==10 | LUCC_map==11 | LUCC_map==12 | LUCC_map==20)=2; 
         LUCC_map(LUCC_map==51 | LUCC_map==52)=3;      %EBF
         LUCC_map(LUCC_map==61 | LUCC_map==62)=4;   %DBF
         LUCC_map(LUCC_map==71 | LUCC_map==72)=5;   %ENF
         LUCC_map(LUCC_map==81 | LUCC_map==82)=6;   %DNF
         LUCC_map(LUCC_map==91 | LUCC_map==92)=7;   %MF
         LUCC_map(LUCC_map==120 | LUCC_map==121 | LUCC_map==122)=8;   %SHR
         LUCC_map(LUCC_map==130)=9;   %GRA
         LUCC_map(LUCC_map==180)=10;   %WET     
         
outputLUCC=LUCC_map;
 end
 
 %% the following functions are for the extracting the LUE data
 function [outputT_timeseries] = get_T_region_openDATA(inputT,location_LAT,location_LON)
              
               LATcentral=(location_LAT(1,1)+location_LAT(1,2))/2;
               LONcentral=(location_LON(1,1)+location_LON(1,2))/2;
               IMAGEstartLON=-180;
               IMAGEstartLAT=90;
               posLON = ceil((LONcentral-IMAGEstartLON)/0.25);
               posLAT = ceil((IMAGEstartLAT- LATcentral)/0.25);
               outputT_timeseries=zeros(60,9);

               
                for num=1:60
                   outputT_timeseries(num,1)=inputT(posLAT-1,posLON-1,num);
                   outputT_timeseries(num,2)=inputT(posLAT-1,posLON,num);
                   outputT_timeseries(num,3)=inputT(posLAT-1,posLON+1,num);
                   outputT_timeseries(num,4)=inputT(posLAT,posLON-1,num);                   
                   outputT_timeseries(num,5)=inputT(posLAT,posLON,num);
                   outputT_timeseries(num,6)=inputT(posLAT,posLON+1,num);                   
                   outputT_timeseries(num,7)=inputT(posLAT+1,posLON-1,num);
                   outputT_timeseries(num,8)=inputT(posLAT+1,posLON,num);
                   outputT_timeseries(num,9)=inputT(posLAT+1,posLON+1,num);
                end
                clear tempFILE

end

function [outputVPD_timeseries] = get_VPD_region_openDATA(inputVPD,location_LAT,location_LON)
               
               LATcentral=(location_LAT(1,1)+location_LAT(1,2))/2;
               LONcentral=(location_LON(1,1)+location_LON(1,2))/2;
               IMAGEstartLON=-180;
               IMAGEstartLAT=90;
               posLON = ceil((LONcentral-IMAGEstartLON)/0.25);
               posLAT = ceil((IMAGEstartLAT- LATcentral)/0.25);
               outputVPD_timeseries=zeros(60,9);
               

               for num=1:60
                   outputVPD_timeseries(num,1)=inputVPD(posLAT-1,posLON-1,num);
                   outputVPD_timeseries(num,2)=inputVPD(posLAT-1,posLON,num);
                   outputVPD_timeseries(num,3)=inputVPD(posLAT-1,posLON+1,num);
                   outputVPD_timeseries(num,4)=inputVPD(posLAT,posLON-1,num);                   
                   outputVPD_timeseries(num,5)=inputVPD(posLAT,posLON,num);
                   outputVPD_timeseries(num,6)=inputVPD(posLAT,posLON+1,num);                   
                   outputVPD_timeseries(num,7)=inputVPD(posLAT+1,posLON-1,num);
                   outputVPD_timeseries(num,8)=inputVPD(posLAT+1,posLON,num);
                   outputVPD_timeseries(num,9)=inputVPD(posLAT+1,posLON+1,num);
               end
                 clear tempFILE
end

function [outputPA_timeseries] = get_PA_region_openDATA(inputPA,location_LAT,location_LON)
%360*540

               LATcentral=(location_LAT(1,1)+location_LAT(1,2))/2;
               LONcentral=(location_LON(1,1)+location_LON(1,2))/2;
               IMAGEstartLON=-180;
               IMAGEstartLAT=90;
               posLON = ceil((LONcentral-IMAGEstartLON)/0.25);
               posLAT = ceil((IMAGEstartLAT- LATcentral)/0.25);
               
               outputPA_timeseries=zeros(60,9);
               
               for num=1:60
                   outputPA_timeseries(num,1)=inputPA(posLAT-1,posLON-1,num);
                   outputPA_timeseries(num,2)=inputPA(posLAT-1,posLON,num);
                   outputPA_timeseries(num,3)=inputPA(posLAT-1,posLON+1,num);
                   outputPA_timeseries(num,4)=inputPA(posLAT,posLON-1,num);                   
                   outputPA_timeseries(num,5)=inputPA(posLAT,posLON,num);
                   outputPA_timeseries(num,6)=inputPA(posLAT,posLON+1,num);                   
                   outputPA_timeseries(num,7)=inputPA(posLAT+1,posLON-1,num);
                   outputPA_timeseries(num,8)=inputPA(posLAT+1,posLON,num);
                   outputPA_timeseries(num,9)=inputPA(posLAT+1,posLON+1,num);
               end
                 clear tempFILE
end


function [outputPAR_timeseries] = get_PAR_region_openDATA(inputPAR,location_LAT,location_LON)

%360*576

               LATcentral=(location_LAT(1,1)+location_LAT(1,2))/2;
               LONcentral=(location_LON(1,1)+location_LON(1,2))/2;
               IMAGEstartLON=-180;
               IMAGEstartLAT=90;
               posLON = ceil((LONcentral-IMAGEstartLON)/0.25);
               posLAT = ceil((IMAGEstartLAT- LATcentral)/0.25);
               
               outputPAR_timeseries=zeros(60,9);
               
                %tempFILE=imread(inputPATH_PAR);
               for num=1:60
                   outputPAR_timeseries(num,1)=inputPAR(posLAT-1,posLON-1,num);
                   outputPAR_timeseries(num,2)=inputPAR(posLAT-1,posLON,num);
                   outputPAR_timeseries(num,3)=inputPAR(posLAT-1,posLON+1,num);
                   outputPAR_timeseries(num,4)=inputPAR(posLAT,posLON-1,num);                   
                   outputPAR_timeseries(num,5)=inputPAR(posLAT,posLON,num);
                   outputPAR_timeseries(num,6)=inputPAR(posLAT,posLON+1,num);                   
                   outputPAR_timeseries(num,7)=inputPAR(posLAT+1,posLON-1,num);
                   outputPAR_timeseries(num,8)=inputPAR(posLAT+1,posLON,num);
                   outputPAR_timeseries(num,9)=inputPAR(posLAT+1,posLON+1,num);
               end
                 clear tempFILE
end
 
 
 
 
 
 
 
 
 %% the following functions are the NDVI rebuilt algorithm


function [rebuiltNDVI,rebuiltQCmat]=rebuiltNDVI_fun(inputNDVI,inputQC)

          effectPIXEL_temp=uint8(inputNDVI>0);
          
          
          effectPIXEL=sum(effectPIXEL_temp,3);
          rebuiltNDVI=uint16(zeros(927,927,60));
          rebuiltQCmat=uint8(zeros(927,927,60));
          
          parfor m=1:927
              for n=1:927
                  if effectPIXEL(m,n)>5
                     [outputNDVIcol,outputQCcol]=gapfilled_NDVI_algorithm(reshape(inputNDVI(m,n,1:60),60,1),reshape(inputQC(m,n,1:60),60,1));
                     rebuiltNDVI(m,n,:)=uint16(outputNDVIcol);
                     rebuiltQCmat(m,n,:)=uint8(outputQCcol);
                     
                  else
                     rebuiltNDVI(m,n,:)=uint16(zeros(60,1));
                     rebuiltQCmat(m,n,:)=uint8(ones(60,1)*3);
                     
                  end
                  %disp(num2str(m));
              end
          end
          %delete(gcp('nocreate'));
end

%  gapfilled by SG filter, maybe is a slow algorithm
function [outputNDVIcol,outputQCcol]=gapfilled_NDVI_algorithm(inputNDVIcol,inputQCcol)
         %inputNDVIcol(1:60,1)=inputNDVIcol_input;
         %outputQCcol(1:60,1)=inputQCcol_input;

        outputQCcol= uint8(inputQCcol);
         
         if inputNDVIcol(1,1)<0
            inputNDVIcol(1,1)=min(inputNDVIcol(inputNDVIcol>0));
         end
         if inputNDVIcol(60,1)<0
            inputNDVIcol(60,1)=min(inputNDVIcol(inputNDVIcol>0));
         end
         
         tags=find(inputNDVIcol>0);
         
         longx=length(tags);         
         outputNDVIcol_filled=inputNDVIcol;

         
         for i=1:longx-1
             tagSTART=tags(i);
             tagEND=tags(i+1);
             % so we have three conditions here 
             %case1 do not need to gapfill, with less than every month have their data, so the tagEND-tagSTART=1
             %case2 need to gapfill but with high quality gapfilled, the data gap less than 4 month, tagEND-tagSTART maxmial is 4 with QC tag ==1  
             %case3 need to gapfill but the result have high uncertainties, we need to have new process method in the future. with QC tag ==2  
             
             if tagEND-tagSTART>=2 && tagEND-tagSTART<=4
                
                for tagNEW=tagSTART+1:tagEND-1
                    outputNDVIcol_filled(tagNEW)=inputNDVIcol(tagSTART)+(inputNDVIcol(tagEND)-inputNDVIcol(tagSTART))*((tagNEW-tagSTART)/(tagEND-tagSTART));
                    outputQCcol(tagNEW)=1;
                end
             elseif tagEND-tagSTART>=4 
                for tagNEW=tagSTART+1:tagEND-1
                    outputNDVIcol_filled(tagNEW)=inputNDVIcol(tagSTART)+(inputNDVIcol(tagEND)-inputNDVIcol(tagSTART))*((tagNEW-tagSTART)/(tagEND-tagSTART));
                    outputQCcol(tagNEW)=2;
                end                
             end
         end


         outputNDVIcol = uint16(sgolayfilt(double(outputNDVIcol_filled),3,5));
         % outputNDVIcol = uint16(outputNDVIcol_filled);
end


 
function [saveTAG]=findlistTAG(filelist)
              filename_temp=filelist(3).name;
              filename_general=filename_temp(1:20);
              counter=0;
              for year=2016:2020
                  for month=1:12
                      counter=counter+1;
                      if month<10
                           filename_temp=strcat(filename_general,num2str(year),'0',num2str(month));
                           
                           [outputTAG]=findTAG(filelist,filename_temp);
                           saveTAG(counter,1)=outputTAG;
                      else
                           filename_temp=strcat(filename_general,num2str(year),num2str(month));
                           [outputTAG]=findTAG(filelist,filename_temp);
                           saveTAG(counter,1)=outputTAG;
                      end
                  end
              end
end


function outputTAG=findTAG(filelist,textname)

         counter=1;
         while(strcmp(filelist(counter).name,textname))==0  && counter<length(filelist)
               counter=counter+1;
         end
         
         
         if counter<length(filelist)
            outputTAG=counter;
         elseif counter==length(filelist)
             if  (strcmp(filelist(counter).name,textname))==1
                 outputTAG=counter;
             else
                 outputTAG=0;
             end
         end

end
