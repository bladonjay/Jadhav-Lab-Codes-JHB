function[out] = setfilterkeepoutput(f, output_name)
% function f = setfilterkeepoutput
% Moves current output to newNameoutput and clears output

for an = 1:length(f)
    f(an).(output_name) = f(an).output;
    f(an).output = [];
end
out = f;


%for i = 1:length(f)
%    tmp_f = f(i);
%    tmp_f = rmfield(f,(outputName));
%    f(i) = tmp_f;
%end
%out = f;
