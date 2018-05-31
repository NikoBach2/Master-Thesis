function progresscomment(comment)
%PROGRESSCOMMENT Show a timestamp and a comment

t = datetime;
timestring = datestr(t, 'HH:MM:SS');
fprintf('%s %s\n', timestring, comment);

end

