% Copyright (c) 2024 Johannes Sieberer
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%%
clear all
%% Get patient files, intial results
listing = dir ("LandMarkReviewer2/*.csv");
cases = struct2table(listing).name;
results.TTTG3D = zeros(height(cases),1);
results = struct2table(results);
m = 0;
%% Iterate through dataset files
for n = 1:height(results)
 %% Import data
    casen = cases(n);
    data = importfile(strcat("(DatsetFolder)/",casen));
    
    %% Get Points (x,y,y)
    femPostCondLat = table2array(data(5,10:12));        %Femoral lateral posterior condyle
    femPostCondMed = table2array(data(6,10:12));        %Femoral medial posterior condyle    
    epicondyleMed = table2array(data(3,10:12));         %Femoral medial epicondyle
    epicondyleLat = table2array(data(4,10:12));         %Femoral lateral epicondyle
    troch = table2array(data(22,10:12));                %Femoral trochlea groove

    tibLatPostCon =  table2array(data(17,10:12));       %Tibia lateral posterior condyle
    tibMedPostCon = table2array(data(18,10:12)) ;       %Tibia medial posterior condyle
    tibMedIterCondTub = table2array(data(16,10:12));    %Tibia medial intracondylar tubercle
    tibLatIterCondTub = table2array(data(15,10:12));    %Tibia lateral intracondylar tubercle
    tibTuber = table2array(data(21,10:12));             %Tibial tuberosity
    tibShaftCen =  table2array(data(19,10:12));         %Tibial shaft center    
    tibBorLat = table2array(data(13,10:12));            %Tibial lateral border
    tibBorMed = table2array(data(14,10:12));            %Tibial medial border

    patLatBod = table2array(data(9,10:12));             %Patella lateral border
    patMedBod = table2array(data(10,10:12));            %Patella medial border
    plateau = table2array(data(20,10:12));              %Tibial plateau
    patDistPol = table2array(data(11,10:12));           %Patella distal pole
    patProxPol = table2array(data(12,10:12));           %Patella proximal pole    


    %% Calculate reference lines
    % Get line shaftcenter to intercondylar midpoint
    tibMidPla = (tibMedIterCondTub + tibLatIterCondTub)/2; %Midpoint between tibial intracondylar tubercles
    tibShaftMid = tibMidPla - tibShaftCen;                 %Tibial long axis (shaft center to midpoint tibial intracondylar tubercles)

    %Get posterior condyle line femur
    femPostCondLine = femPostCondLat -femPostCondMed;      %Line between femoral posterior condyles
    femCondLineProj = femPostCondLine - dot(femPostCondLine,tibShaftMid)*tibShaftMid/norm(tibShaftMid)^2; %Femoral condyle line perpendicular to the tibial longitudinal axis
    femCondLineProjNorm = femCondLineProj / norm(femCondLineProj); %Direction vector of the femoral post condyle line
    
    %Get posterior condyle line tibia
    tibPostCondLine = tibLatPostCon - tibMedPostCon;       %Line between tibial posterior condyles
    tibCondLineProj = tibPostCondLine - dot(tibPostCondLine,tibShaftMid)*tibShaftMid/norm(tibShaftMid)^2;% Tibial condyle Line perpendicular to the tibial longitudinal axis
    tibCondLineProjNorm = tibCondLineProj / norm(tibCondLineProj); %Direction of tibial post condyle line
    
    %% Tibiofemoral rotation in degree
    rottf= signedAngleTwo3DVectors(femCondLineProj,tibCondLineProj,tibShaftMid,1);%Angle between posterior condyle lines, along tibial long axis
    results.ROTTF(n) = rottf*sign(femCondLineProjNorm(1))*180 / pi;               %conversion to degree and sign for external and interal rotation (left/right knee)
    
    %% 3D TT-TG 
    tttg = tibTuber - troch;   %Line between trochlea groove and tibial tuberosity
    results.TTTG3D(n) = dot(femCondLineProjNorm, tttg)/norm(femCondLineProjNorm); %Projection on the posterior condyle line
   
    %% 2D TT-TG
    femPostCondLineYX = [femPostCondLine(1), femPostCondLine(2), 0]; %Projection of femoral posterior condyle line onto axial plane
    results.TTTG2D(n) = dot(femPostCondLineYX, tttg)/norm(femPostCondLineYX); %Projection of TT-TG line on axial posterior condyle line (classic 2D TTTG)
    
    %% Patellar tilt
    patLineXY = patMedBod(1:2) - patLatBod(1:2); %%Get axial projection of the patella line from lateral to medial border

    if(femPostCondLineYX(1) > 0) %% Get angle sign right and claculate angle between axial femoral posterior condyle line and patella line
        results.Patellartilt(n) = signedAngleTwo3DVectors([-femPostCondLineYX],[patLineXY,0],[0,0,1],1)*180/pi;
    else
        results.Patellartilt(n) = signedAngleTwo3DVectors([-femPostCondLineYX],[patLineXY,0],[0,0,1],1)*180/pi * -1; 
    end
    % Signed two angles credits: Seth Wagenman (2024). signedAngleTwoVectors (https://www.mathworks.com/matlabcentral/fileexchange/78300-signedangletwovectors), MATLAB Central File Exchange. Retrieved November 5, 2024. 

    %% Patellar length
    results.PatLength(n) = norm(patDistPol - patProxPol);

    %% Patellar height
    results.PatHeight(n) = norm(patDistPol-plateau) / results.PatLength(n);  %%Ratio between distance between the tibial plateau and the distal pat pole to patella length
    
    %% Epricondylar distance
    epi = epicondyleMed - epicondyleLat; %Line between femoral epicondyles
    results.FemEpiCond(n) =dot(femCondLineProjNorm, epi); %Length of epicondyle line projected on the posterior condyle line

    %% Tibial plateau width
    tpw = tibBorLat - tibBorMed; %Line between tibial pleateau borders
    results.TibPlatWidth(n) = dot(tibCondLineProjNorm, tpw); %Length of tibial border line projected on the posterior condyle line

    %% TT-TIM Tibial tuberosity to tibial intercondyal midpoint
    %Calcualte mid point per orignal paper figure (midpoint of bounding box tibial lateral and medial borders and tibial posterior condyles plus anterior tibial plateau 
    tibLatConPl =  tibLatPostCon - plateau;
    tibplaMidVect = (tibLatConPl - dot(tibCondLineProjNorm,tibLatConPl)*tibCondLineProjNorm)/2;
    tpwMidVect = dot(tibCondLineProjNorm, tpw)*tibCondLineProjNorm/2;
    
    tibLatConLatBor = tibBorLat - tibLatPostCon;
    tibLatConLatBorProjVec = tibLatConLatBor - dot(tibCondLineProjNorm,tibLatConLatBor)*tibCondLineProjNorm;

    midPointPl = tibBorLat - tibLatConLatBorProjVec + tibplaMidVect - tpwMidVect; %Midpoint plateau

    tttm = tibTuber - midPointPl; %Tibial tuberosity to mid piont tibial plateau
    results.TTTIM(n) = dot(tibCondLineProjNorm,tttm); %Distance projected onto the tibial posterior condyle line
    
    %% TT-TEM Tibial tuberosity to tibial eminence midpoint
    tttem = tibTuber - tibMidPla; %%Line between TT and tib plateau midpoint
    results.TTTEM(n) = dot(tibCondLineProjNorm, tttem)/norm(tibCondLineProjNorm); %%Projection on tibial posterior condyle line

    %% Sagittal TT-TG
    sttg = tttg - dot(femCondLineProjNorm, tttg)*femCondLineProjNorm; %%Projection on sagittal plane
    results.sTTTG(n) = norm(sttg(1:2))*sign(sttg(2)); %sign for posterior and anterior to trochlea(+/-)
end
