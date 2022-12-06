function animalinfo = animaldef(animalname)

switch animalname
    
    % Ripple Disruption Expt Animals
        
    case 'sjc'
        animalinfo = {'sjc', '/data25/sjadhav/RippleInterruption/sjc_direct/', 'sjc'};
    case 'RE1'
        animalinfo = {'RE1', '/data25/sjadhav/RippleInterruption/RE1_direct/', 'RE1'};
    case 'RNa'
        animalinfo = {'RNa', '/data25/sjadhav/RippleInterruption/RNa_direct/', 'RNa'};
    case 'RNb'
        animalinfo = {'RNb', '/data25/sjadhav/RippleInterruption/RNb_direct/', 'RNb'};
    case 'RNc'
        animalinfo = {'RNc', '/data25/sjadhav/RippleInterruption/RNc_direct/', 'RNc'};
    case 'RNd'
        animalinfo = {'RNd', '/data25/sjadhav/RippleInterruption/RNd_direct/', 'RNd'};
    case 'RCa'
        animalinfo = {'RCa', '/data25/sjadhav/RippleInterruption/RCa_direct/', 'RCa'};
    case 'RCb'
        animalinfo = {'RCb', '/data25/sjadhav/RippleInterruption/RCb_direct/', 'RCb'};
    case 'RCc'
        animalinfo = {'RCc', '/data25/sjadhav/RippleInterruption/RCc_direct/', 'RCc'};
    case 'RCd'
        animalinfo = {'RCd', '/data25/sjadhav/RippleInterruption/RCd_direct/', 'RCd'};
    case 'REc'
        animalinfo = {'REc', '/data25/sjadhav/RippleInterruption/REc_direct/', 'REc'};
    case 'REd'
        animalinfo = {'REd', '/data25/sjadhav/RippleInterruption/REd_direct/', 'REd'};
    case 'REe'
        animalinfo = {'REe', '/data25/sjadhav/RippleInterruption/REe_direct/', 'REe'};
    case 'REf'
        animalinfo = {'REf', '/data25/sjadhav/RippleInterruption/REf_direct/', 'REf'};
    case 'REg'
        animalinfo = {'REg', '/data25/sjadhav/RippleInterruption/REg_direct/', 'REg'};
    case 'REh'
        animalinfo = {'REh', '/data25/sjadhav/RippleInterruption/REh_direct/', 'REh'};
    case 'HPa'
        animalinfo = {'HPa', '/data25/sjadhav/HPExpt/HPa_direct/', 'HPa'};
    case 'HPb'
        animalinfo = {'HPb', '/data25/sjadhav/HPExpt/HPb_direct/', 'HPb'};
    case 'HPc'
        animalinfo = {'HPc', '/data25/sjadhav/HPExpt/HPc_direct/', 'HPc'};
    case 'Ndl'
        animalinfo = {'Ndl', '/data25/sjadhav/HPExpt/Ndl_direct/', 'Ndl'};
    case 'Nadal'
        animalinfo = {'Ndl', '/data25/sjadhav/HPExpt/Ndl_direct/', 'Ndl'};
    case 'Rtl'
        animalinfo = {'Rtl', '/data25/sjadhav/HPExpt/Rtl_direct/', 'Rtl'};
    case 'Rosenthal'
        animalinfo = {'Rtl', '/data25/sjadhav/HPExpt/Rtl_direct/', 'Rtl'};
    case 'Brg'
        animalinfo = {'Brg', '/data25/sjadhav/HPExpt/Brg_direct/', 'Brg'};
    case 'Borg'
        animalinfo = {'Brg', '/data25/sjadhav/HPExpt/Brg_direct/', 'Brg'};
    case 'CS31'
        animalinfo = {'CS31','E:\Brandeis datasets\OdorPlaceAssociation\CS31Expt\CS31_direct\','CS31'};
    case 'CS32'
        animalinfo = {'CS32','E:\Brandeis datasets\OdorPlaceAssociation\CS32Expt\CS32_direct\','CS32'};
        
    case 'CS33'
        animalinfo = {'CS33','E:\Brandeis datasets\OdorPlaceAssociation\CS33Expt\CS33_direct\','CS33'};
    case 'CS34'
        animalinfo = {'CS34','E:\Brandeis datasets\OdorPlaceAssociation\CS34Expt\CS34_direct\','CS34'};
    case 'CS35'
        animalinfo = {'CS35','E:\Brandeis datasets\OdorPlaceAssociation\CS35Expt\CS35_direct\','CS35'};
    case 'CS36'
        animalinfo = {'CS36','E:\Brandeis datasets\OdorPlaceAssociation\CS36Expt\CS36_direct\','CS36'};
    case 'CS39'
        animalinfo = {'CS39','E:\Brandeis datasets\OdorPlaceAssociation\CS39Expt\CS39_direct\','CS39'};
    case 'CS41'
        animalinfo = {'CS41','E:\Brandeis datasets\OdorPlaceAssociation\CS41Expt\CS41_direct\','CS41'};
    case 'CS42'
        animalinfo = {'CS42','E:\Brandeis datasets\OdorPlaceAssociation\CS42Expt\CS42_direct\','CS42'};
    case 'CS44'
        animalinfo = {'CS44','E:\Brandeis datasets\OdorPlaceAssociation\CS44Expt\CS44_direct\','CS44'};
        
    otherwise
        
        error(['Animal ',animalname, ' not defined.']);
end
