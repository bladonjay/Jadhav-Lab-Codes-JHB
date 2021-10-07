function [opts] = setStateScriptOpts()

 opts = delimitedTextImportOptions("NumVariables", 2);
    
    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = "#";
    
    % Specify column names and types
    opts.VariableNames = ["VarName1", "Var2"];
    opts.SelectedVariableNames = "VarName1";
    opts.VariableTypes = ["string", "string"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Specify variable properties
    opts = setvaropts(opts, ["VarName1", "Var2"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["VarName1", "Var2"], "EmptyFieldRule", "auto");
end

