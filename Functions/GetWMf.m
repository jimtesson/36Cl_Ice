function [ ] = GetWMf( name )
% Function that first applies the Chauvenet criterion to a set of data and
% then compute a MSWD-weighted mean out of filtered dataset

% Open paths and load data
addpath Functions
addpath Constants
addpath jsonlab

% Retreive json file
Data=loadjson(name);
% Name=fieldnames(Data);
% Data=getfield(Data,Name{1});

% Load values
PRval=Data.PRval;
PRerr=Data.PRerr;
NbSig=Data.NbSig;
Chauv=Data.Chauv;
Wstdval=Data.Wstdval;

% Calculate
% NbSig=2;
if Chauv ==1;
    [PRval,PRerr]=Chauvenet(PRval,PRerr,NbSig);
end
[~,Wmean,ErrWm]=CorrWM([PRval PRerr]);

% Different error bar
if Wstdval==1 ;
	if length(PRval) > 1;
   		ErrWm = Wstd (PRval, (1./PRerr).^2);
	end
end

% Export
DataOut.WM=[Wmean ErrWm];

% Write json
DataOut=savejson(name,DataOut);
NameOut=strcat(name(1:end-2),'out');
fileID=fopen(NameOut,'w');
fprintf(fileID,'%s',DataOut);
fclose(fileID);

end

