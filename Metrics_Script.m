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
casen = cases(n-floor((m+1)/2));
data = importfile(strcat("LandMarkReviewer2/",casen));
    if(height(data) <= 23)
        %% Get Points (x,y,y)
        tibLatPostCon =  table2array(data(16,10:12));       %tibia lateral posterior condyle
        tibMedPostCon = table2array(data(17,10:12)) ;       %tibia medial  posterior condyle
        tibMedIterCondTub = table2array(data(15,10:12));    %tibia lateral intercondylar tubercle
        tibLatIterCondTub = table2array(data(14,10:12));    %tibia medial  intercondylar tubercle
        tibTuber = table2array(data(20,10:12));             %tibial tuberosity
        tibShaftCen =  table2array(data(18,10:12));         %tibial shaft center
        femPostCondLat = table2array(data(4,10:12));        %femoral lateral posterior condyle
        femPostCondMed = table2array(data(5,10:12));        %femoral medial  posterior condyle
        patLatBod = table2array(data(8,10:12));             %patella lateral border
        patMedBod = table2array(data(9,10:12));            %patella medial  border
        troch = table2array(data(21,10:12));                %trochlea groove
        plateau = table2array(data(19,10:12));              %tibial plateau
        patDistPol = table2array(data(10,10:12));           %patella distal   pole
        patProxPol = table2array(data(11,10:12));           %patella proximal pole
        results.Cases(n) = casen; %assign case number
    end
    
    %% Calculate common vectors
    % Get line tib shaftcenter to intercondylar midpoint
    tibMidPla = (tibMedIterCondTub + tibLatIterCondTub)/2;
    tibShaftMid = tibMidPla - tibShaftCen; %Tibial longitudinal axis
    tibShaftNorm = tibShaftMid/norm(tibShaftMid)^2;
    % Get posterior  femur and tibial condyle lines
    femPostCondLine = femPostCondLat -femPostCondMed;
    femCondLineProj = femPostCondLine - dot(femPostCondLine,tibShaftNorm)*tibShaftNorm; 
    femCondLineProjNorm  = femCondLineProj/norm(femCondLineProj); %Projected femoral condyle line normalized
    tibPostCondLine = tibLatPostCon - tibMedPostCon;
    tibCondLineProj = tibPostCondLine - dot(tibPostCondLine,tibShaftNorm)*tibShaftNorm; 
    tibCondLineProjNorm =  tibCondLineProj/norm(tibCondLineProj); %Projected tibial condyle line normalized
    % Patella line through medial and posterior poles
    patPoleLine = patMedBod - patLatBod;
    %Tuberosity to tibial plateau eminence
    tttp =  tibMidPla - tibTuber;
    %Trochlea groove to tibial plateau eminence
    tgtp =  troch - tibMidPla;
    %Trochlea groove to tibial tuberosity
    tttg = tibTuber - troch;
    
    %% Calculate Tibiofemoral rotation
    rotVector = cross(femCondLineProj,tibCondLineProj); %Assume knee is generally superior inferior aligned (third component pointing up and down)
    results.ROTTF(n) = signedAngleTwo3DVectors(femCondLineProj,tibCondLineProj,rotVector(:) .* sign(rotVector(3)).* sign(femCondLineProj(1)),1)*180 / pi;

    %% 3D TT-TG
    results.TTTG3D(n) = dot(femCondLineProjNorm, tttg);
    
    %% 2D TT-TG
    femPostCondLineYX = [femPostCondLine(1), femPostCondLine(2), 0]; %Femur posterior condyle line projected on axial slice
    XYTTTG = [tibTuber(1) - troch(1),tibTuber(2) - troch(2),0];
    results.TTTG2D(n) = dot(femPostCondLineYX, XYTTTG)/norm(femPostCondLineYX);
    
    %% Patella tilt
    femLineXY = -[femPostCondLine(1),femPostCondLine(2)];    %Axial slice projection femoral condyle line
    patLineXY = [patPoleLine(1),patPoleLine(2)];            %Axial slice projection patella line
    if(femLineXY(1) < 0) %Left and right knees
    results.Patellartilt(n) = signedAngleTwo3DVectors([femLineXY,0],[patLineXY,0],[0,0,1],1)*180/pi;
    else
    results.Patellartilt(n) = signedAngleTwo3DVectors([femLineXY,0],[patLineXY,0],[0,0,1],1)*180/pi * -1; 
    end
    
    %% Tuberosity distance    
    results.TubDist(n) = abs(dot(tibCondLineProjNorm, tttp)/norm(tibCondLineProjNorm));
    
    %% Trochlea groove distance
    results.TGDist(n) = -dot(femCondLineProjNorm, tgtp)/norm(femCondLineProjNorm);

    %% Translational TT-TG
    results.TransTTTG(n) =  results.TGDist(n) + results.TubDist(n);

    %% Rotational TT-TG
    results.RotDist(n) = results.TTTG3D(n) - results.TransTTTG(n) ;

    %% Patella Height
    results.PatHeight(n) = norm(patDistPol-plateau) / norm(patProxPol - patDistPol);
    
    %% Sagittal TT-TG
    sttg = tttg - dot(femCondLineProjNorm, tttg)*femCondLineProjNorm - dot(tttg,tibShaftNorm)*tibShaftNorm;
    results.sTTTG(n) = norm(sttg)*sign(sttg(2));
end
