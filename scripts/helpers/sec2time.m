function[time] = sec2time(sec)
    time = datestr(toc/(24*60*60),'HH:MM:SS');
end