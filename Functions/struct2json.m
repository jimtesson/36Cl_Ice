function [ ] = struct2json( struct, outname )
% Create a json from the struct structure

struct=mat2json(struct);
fileID = fopen(outname,'w');
fprintf(fileID,'%s',struct);
fclose(fileID);

end

