function StrippedString = stripz(number)
    %Uses a regular expression search to strip lagging zeroes from a number
    %which is put into fixed-point format
    
    str = sprintf('%.14f ',number);
    if floor(number)==number
        str = regexprep(str, '\.[0]+ ', '');
    else
        str = regexprep(str, '[0]+ ', '');
    end
    StrippedString = str;

end