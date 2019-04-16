
% command eeglab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; 
home_path = [pwd() '/../'];
subject_list = {'pilot1', 'pilot2', 'pilot3', 'pilot4'};
%subject_list = {'pilot1', 'pilot2'};

% flag for selecting the window and step for artifcat rejection

window = 0
step = 1

% open a file to save erpset_list
f = fopen('erpset_list.txt','w');

%% Data Analysis 
for s=1:length(subject_list)

    %% EEGLAB operation
    % setting a variable to use often
    data_path = [home_path 'Data/' subject_list{s} '/'];
    
    % load brain vision data
    EEG = pop_loadbv([data_path 'rawdata/'], [subject_list{s} '.vhdr']);
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

    %% Start of the second part
    % epoch
    EEG = pop_epochbin( EEG , [-500.0  500.0],  'pre'); 
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
                              50 );
        EEG.setname = [EEG.setname '_arst'];
        EEG = pop_saveset(EEG, 'filename',[EEG.setname '.set'],'filepath', data_path);
    end;
    
    % clear artifact
    EEG  = pop_resetrej( EEG , 'ResetArtifactFields', 'on' ); 
    EEG.setname = [EEG.setname '_resetrej'];
    EEG = pop_saveset(EEG, 'filename',[EEG.setname '.set'],'filepath', data_path);

    EEG_faces = EEG;
    EEG_cars = EEG;
    
    % Create epoched dataset for the first condition (faces)
    % EEG = pop_loadset('filename',[subject_list{s} filname_ext], 'filepath', data_path);
    EEG_faces = pop_epoch( EEG_faces, {  'B1(S1)'  'B1(S10)'  'B1(S11)'  'B1(S12)'  'B1(S13)'  'B1(S14)'  'B1(S15)'  'B1(S16)'  'B1(S17)'  'B1(S18)'...
                            'B1(S19)'  'B1(S2)'  'B1(S20)'  'B1(S21)'  'B1(S22)'  'B1(S23)'  'B1(S24)'  'B1(S25)'  'B1(S26)'  'B1(S27)'...
                            'B1(S28)'  'B1(S29)'  'B1(S3)'  'B1(S30)'  'B1(S31)'  'B1(S32)'  'B1(S33)'  'B1(S34)'  'B1(S35)'  'B1(S36)'...
                            'B1(S37)'  'B1(S38)'  'B1(S39)'  'B1(S4)'  'B1(S40)'  'B1(S5)'  'B1(S6)'  'B1(S7)'  'B1(S8)'  'B1(S9)'   }, [-0.5 0.5]);
    EEG_faces = pop_rmbase( EEG_faces, [-500    0]);
    EEG_faces.setname = [subject_list{s} '_bin1_faces'];
    EEG_faces = pop_saveset(EEG_faces, 'filename', [EEG_faces.setname '.set'], 'filepath', data_path);

    % Create epoched dataset for the first condition (cars) I need to change bins.     
    EEG_cars = pop_epoch( EEG_cars, {  'B1(S1)'  'B1(S10)'  'B1(S11)'  'B1(S12)'  'B1(S13)'  'B1(S14)'  'B1(S15)'  'B1(S16)'  'B1(S17)'  'B1(S18)'...
                            'B1(S19)'  'B1(S2)'  'B1(S20)'  'B1(S21)'  'B1(S22)'  'B1(S23)'  'B1(S24)'  'B1(S25)'  'B1(S26)'  'B1(S27)'...
                            'B1(S28)'  'B1(S29)'  'B1(S3)'  'B1(S30)'  'B1(S31)'  'B1(S32)'  'B1(S33)'  'B1(S34)'  'B1(S35)'  'B1(S36)'...
                            'B1(S37)'  'B1(S38)'  'B1(S39)'  'B1(S4)'  'B1(S40)'  'B1(S5)'  'B1(S6)'  'B1(S7)'  'B1(S8)'  'B1(S9)'   }, [-0.5 0.5]);
    EEG_cars = pop_rmbase( EEG_cars, [-500    0]);
    EEG_cars.setname = [subject_list{s} '_bin2_cars'];
    EEG_cars = pop_saveset(EEG_cars, 'filename', [EEG_cars.setname '.set'], 'filepath', data_path);

    
    % save erpset name onto the text file "erpset_list.txt"
    fprintf(f,[data_path subject_list{s} '_diff_ERPs.erp\n']);

end;

% close the file
fclose(f);



