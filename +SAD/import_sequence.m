function sequence = import_sequence(folder)

content = dir(folder);
if length(content) > 0
	i = 1;
	while i < length(content) + 1
   	 isfasta = strfind(content(i).name,'.fasta');
	    if isfasta
   	     break;
	    end
   	 i = i+1;
	end
	if ~isfasta  
   	 fprintf('Error: Fasta sequence file not found in %s folder',folder);
	end
	initread = fastaread([folder content(i).name]);
	sequence = initread.Sequence;
else
	fprintf('The folder %s does not exist\n.',folder);
end
