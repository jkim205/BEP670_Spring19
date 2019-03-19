
% command eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; 
home_path = [pwd() '/../'];
% subject_list = {'pilot1', 'pilot2', 'pilot3', 'pilot4'};
subject_list = {'pilot1'};

window = 1
step = 0

%% Data Analysis 
for s=1:length(subject_list)

    % setting a variable to use often
    data_path = [home_path 'Data/' subject_list{s} '/']
    
    % load brain vision data
    EEG = pop_loadbv([data_path '/rawdata/'], [subject_list{s} '.vhdr']);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'gui','off'); 
    EEG = pop_editset(EEG, 'setname', subject_list{s});
    EEG = pop_saveset(EEG, 'filename',[EEG.setname '.set'],'filepath', data_path);
    
    % plot eeg file
    % pop_eegplot(EEG, 1, 1, 1);
    
    % run basic filter
    EEG  = pop_basicfilter( EEG,  1:31 , 'Boundary', 'boundary', 'Cutoff', [ 0.1 30], 'Design', 'butter', 'Filter', 'bandpass', 'Order',2 );
    EEG.setname = [EEG.setname, '_filt'];
    EEG = pop_saveset(EEG, 'filename',[EEG.setname '.set'],'filepath', data_path);

    % plot eeg file
    % pop_eegplot(EEG, 1, 1, 1);
    
    % create eventlist
    EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' },...
           'Eventlist', [data_path 'elist.txt']);
    EEG.setname = [EEG.setname '_elist'];
    EEG = pop_saveset(EEG, 'filename',[EEG.setname '.set'],'filepath', data_path);

    % plot eeg file
    % pop_eegplot(EEG, 1, 1, 1);
    
    % create bins
    EEG  = pop_binlister(EEG , 'BDF', [home_path 'BinFiles/binDescriptor.txt'], 'IndexEL',  1, 'SendEL2', 'EEG', 'Voutput', 'EEG' );
    EEG.setname = [EEG.setname '_bins'];
    EEG = pop_saveset(EEG, 'filename',[EEG.setname '.set'],'filepath', data_path);

    % epoch
    EEG = pop_epochbin( EEG , [-100.0  500.0],  'pre'); 
    EEG.setname = [EEG.setname '_be'];
    EEG = pop_saveset(EEG, 'filename',[EEG.setname '.set'],'filepath', data_path);        
       
    % artifact rejection - moving window to window
    if(window)
        EEG  = pop_artmwppth( EEG , 'Channel',  1:16, 'Flag',  1, 'Threshold',  100, 'Twindow', [ -100 498], 'Windowsize',  200, 'Windowstep',  100 ); 
        EEG.setname = [EEG.setname '_arw'];
        EEG = pop_saveset(EEG, 'filename',[EEG.setname '.set'],'filepath', data_path);
    end;

    % artifact rejection - step
    if(step)
        EEG  = pop_artstep( EEG , 'Channel',  1:16, 'Flag',  1, 'Threshold',  100, 'Twindow', [ -100 498], 'Windowsize',  200, 'Windowstep',...
                              50 )
        EEG.setname = [EEG.setname '_arst'];
        EEG = pop_saveset(EEG, 'filename',[EEG.setname '.set'],'filepath', data_path);
    end;
    
    % clear artifact
    EEG  = pop_resetrej( EEG , 'ResetArtifactFields', 'on' ); 
    EEG.setname = [EEG.setname '_resetrej'];
    EEG = pop_saveset(EEG, 'filename',[EEG.setname '.set'],'filepath', data_path);

    % make new erpset
    ERP = pop_averager( EEG , 'Criterion', 'good', 'ExcludeBoundary', 'on', 'SEM', 'on' );
    ERP = pop_savemyerp( ERP, 'erpname', [subject_list{s} '_ERPs.set'], 'filename', [subject_list{s} '_ERPs.erp'],'filepath', data_path);
    
    % plot
    % pop_ploterps( ERP, 1:4 ,1:31)
    
    % bin operation
    ERP = pop_binoperator( ERP, {  'b5 = b1-b3 label Faces minus sFaces',...
  'b6 = b2-b4 label Cars minus sCars', 'b7 = b5-b6 label DFaces minus DCars', 'b8=b1-b2 label Faces minus Cars'});
    ERP = pop_savemyerp(ERP, 'erpname', [subject_list{s} 'diff_ERPs.set'], 'filename', [subject_list{s} 'diff_ERPs.erp'], 'filepath',...
 data_path, 'Warning', 'on');

    % plot
    pop_ploterps(ERP, 7:8 ,1:31)
    
    % measurement
    ERP = pop_geterpvalues( ERP, [ 140 210], [ 5 6], [ 14 15 19 20] , 'Baseline', 'pre', 'FileFormat', 'wide', 'Filename', [data_path 'measurement.txt'],...
 'Fracreplace', 'NaN', 'InterpFactor',  1, 'Measure', 'meanbl', 'PeakOnset',  1, 'Resolution',  3 );
    
end;

