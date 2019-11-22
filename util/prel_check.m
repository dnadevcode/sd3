function check = prel_check(experiment,functions,actions);

    folderEndsCorr = strcmp(experiment.targetFolder(end),'/');
    if folderEndsCorr
        check = 1;
    else
        check = 0;
        print_errors(folderEndsCorr);
    end
end

function print_errors(folderEndsCorr)
    if ~folderEndsCorr
        fprintf('\n\n INPUT ERROR: Designated "targetFolder" must end in "/".\n\n')
    end
end
