% This function parses the information of the parameter xml-files used for the cpp HJBQVI_solve
% into a matlab array
function s = parseFEniCSxml(filename, p)

try
   tree = xmlread(strcat('xml/',filename));
catch
   error('Failed to read XML file %s.',filename);
end
listOfParameters = tree.getElementsByTagName('parameter');
N = listOfParameters.getLength;
if nargin < 2
	s = struct;
else
	s = p;
end 

for k = 1:N
	item =listOfParameters.item(k-1); % indexing begins with 0 instead of 1
	% keys(k)   = item.getAttribute('key');
	% types(k)  = item.getAttribute('type');
	% values(k) = item.getAttribute('value');
	ename = char(item.getAttribute('key'));
	etype = char(item.getAttribute('type'));
	evalue = char(item.getAttribute('value')); % java.lang.Strings can only be converted to chars directly
	if strcmp(etype,'double')
		evalue = str2double(evalue);
	else
		if strcmp(etype,'bool')
			if strcmp(evalue, 'true')
				evalue = 1;
			else
				evalue = 0;
			end % if strcmp(evalue, 'true')
		else
			if strcmp(etype,'int')
				evalue = str2double(evalue);
			else
				if strcmp(etype,'string')
					evalue = string(evalue);
				else
					error('type of element in xml-file not known')
				end % if strcmp(etype,'string')
			end % if strcmp(etype,'int')
		end % if strcmp(etype,'bool')
	end % if strcmp(etype,'double')

	% Add element to structure array
	s.(ename) = evalue;

end % for k = 1:N

end % function s = parseFEniCSxml(filename)
