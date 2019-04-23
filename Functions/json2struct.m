function [ Struct ] = json2struct( name )
% Turn a json file in the current directory (or connected directories) into
% a matlab structure. 

fileID=fopen(name,'r');
Jsontxt=fscanf(fileID,'%s');
% Data=JSON.parse(Jsontxt); % 2 options possible
Struct=json2mat(Jsontxt);

end

