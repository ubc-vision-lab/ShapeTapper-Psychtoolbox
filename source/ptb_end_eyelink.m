function ptb_end_eyelink(edfFile, edf_file_out)
%ptb_end_eyelink Ends EyeLink recording for this experiment instance

WaitSecs(0.1);
Eyelink('CloseFile');
% download data file
try
    fprintf('Receiving data file ''%s''\n', edfFile );
    status=Eyelink('ReceiveFile');
    if status > 0
        fprintf('ReceiveFile status %d\n', status);
    end
    if 2==exist (edfFile, 'file')
        fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd);
    end
catch
    fprintf('Problem receiving data file ''%s''\n', edfFile );
end
Eyelink('ShutDown');

movefile(edfFile,edf_file_out);

end

