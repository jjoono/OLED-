function C = csv2cell_simple(f)
L = strsplit(strtrim(fileread(f)), {'\n','\r\n'}); C = {};
for i = 2:numel(L)
  if isempty(strtrim(L{i})), continue; end
  C(end+1,:) = strsplit(strtrim(L{i}), ',');
end
end
