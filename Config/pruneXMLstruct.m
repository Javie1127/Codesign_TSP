function S_out = pruneXMLstruct(S_in)
% function to further prune struct loaded using xml2struct
tmpS = S_in;
if (isstruct(S_in))
    childNodeNames = fieldnames(S_in);
    numChildren = numel(childNodeNames);
    for i = 1:numChildren
        if (strcmp(childNodeNames{i}, 'Text'))
            % Remove the "Text" layer
            tmpS = S_in.(genvarname(childNodeNames{i}));
            % Convert to number if parameter is a number
            [X,is_param_num]=str2num(tmpS);
            if is_param_num
                tmpS = X;
            end
            % tmpS
            if (numChildren>1)
                error('XML Config Error: Text can only be defined at the leaf node!\n');
            end
            if (isempty(tmpS))
                error('XML Config Error: Unconfigured Parameter in XML!\n');
            end
        else
            tmpS.(genvarname(childNodeNames{i})) = pruneXMLstruct(S_in.(genvarname(childNodeNames{i})));
        end
    end
end
S_out = tmpS;
end