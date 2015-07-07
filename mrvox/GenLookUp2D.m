function Dico = GenLookUp2D(ModelFileName,SeqFileName,PathOut,ServerInfo)

% This function generate the LookUp table defined in FileIn and save as a file in PathOut.
%
% Input: - FileIn: path to the file where the Model and the MR Sequence are defined
%               (Must be a a cell of strings for multiple files)
%        - PathOut: path to the repertory where the dictionaries will be
%        saved
% Optional: - ServerInfo: Path to the file defining the cluster info for
% parallel jobs (required MCDS + Matlab Runtime Library installed)

if nargin < 4
    Para = 0;
    if nargin == 2
        PathOut = '';
        ServerInfo = '';
    elseif nargin == 3
        ServerInfo = '';
    end
elseif nargin == 4
    Para = 1;
end

if ~iscell(ModelFileName)
    FileIn = {ModelFileName};
    NbFiles = 1;
else
    FileIn = ModelFileName;
    NbFiles = numel(ModelFileName);
end

for NbFile = 1:numel(FileIn)
    clear results Dico JobSeq
    fprintf('Generate Dictionary for File: %d/%d\n', NbFile,NbFiles);
    File = FileIn{NbFile};
    
    %% Read input parameters
    Model   = ReadModel(File);
    Seq     = SelectSeq(SeqFileName,1);
    [Label, par]    = GenerateInputDico(Model);
    
    NbSimu = size(par,1);
    fprintf('Number of elements in the dictionary: %d\n',NbSimu);
    
    %% Run Parallel Jobs on the cluster
    if Para
        Cluster = ReadClusterInfo(ServerInfo);
        
        %% Split in different Jobs to avoid overloading
        if NbSimu > Cluster.MaxTasksPerJob
            remain = mod(NbSimu,Cluster.MaxTasksPerJob);
            for aa=1:((NbSimu-remain)/(Cluster.MaxTasksPerJob))
                JobSeq(aa,:) = [(aa-1)*Cluster.MaxTasksPerJob+1 aa*Cluster.MaxTasksPerJob];
            end
            if remain > 0
                JobSeq(aa+1,:) = [JobSeq(aa,2)+1  JobSeq(aa,2)+remain];
            end
        else
            JobSeq = [1 size(par,1)];
        end
        
        %% Submit Jobs
        % Find scheduler
        jm = findResource('scheduler','type',Cluster.Scheduler.Type,'Name',Cluster.Scheduler.Name);
        if isempty(jm);
            fprintf('Could not find the server\n Exiting\n');
            return
        end
        
        % Submit jobs
        for JobId = 1:size(JobSeq,1)
            job(JobId) = createJob(jm);
            set(job(JobId),'MaximumNumberOfWorkers',Cluster.NbOfWorkers);
            set(job(JobId),'PathDependencies',Cluster.FilePath);
            fprintf('Block Input from %d/%d to %d/%d\n',JobSeq(JobId,1),NbSimu,JobSeq(JobId,2),NbSimu);
            for SimuIdx = JobSeq(JobId,1):JobSeq(JobId,2)
                if Cluster.MRL.On % Used the standalone version with the MRL routine
                    fprintf('Not implemented yet. Exiting\n')
                    return
                else % Used the matlab version
                    createTask(job(JobId), @VoxelSim2D_do_one, 2, {File,SeqFileName,SimuIdx,1});
                end
            end
            submit(job(JobId));
            waitForState(job(JobId))
            results(JobSeq(JobId,1):JobSeq(JobId,2),:) = getAllOutputArguments(job(JobId));
            
            % Check for errors
            task_failed = [];
            for tas = 1:numel(job(JobId).Task)
                if ~isempty(job(JobId).Tasks(tas).ErrorIdentifier)
                    task_failed = [task_failed tas+JobSeq(JobId,1)-1];
                end
            end
            
            % re-run if errors
            if ~isempty(task_failed)
                task_failed_remain = task_failed;
                try_n_times_max = 3;
                try_n_times = try_n_times_max;
                while ~isempty(task_failed_remain) && (try_n_times > 0)
                    fprintf(['Retrying ' sprintf('(%d/%d) ', try_n_times_max - try_n_times+1, try_n_times_max) ...
                        'failed tasks: ' sprintf('%d ', task_failed) '\n']);
                    job_failed = createJob(jm);
                    for tas = 1:numel(task_failed_remain)
                        createTask(job_failed, @VoxelSim2D_do_one, 2, {File, SeqFileName, task_failed_remain(tas), 1});
                    end
                    submit(job_failed);
                    waitForState(job_failed);
                    res_failed = getAllOutputArguments(job_failed);
                    task_failed = task_failed_remain;
                    task_failed_remain = [];
                    for tas = 1:numel(job_failed.Task)
                        if isempty(job_failed.Tasks(tas).ErrorIdentifier) % task succeeded
                            results(task_failed(tas),:) = res_failed(tas, :);
                        else
                            task_failed_remain = [task_failed_remain task_failed(tas)];
                        end
                    end
                    try_n_times = try_n_times - 1;
                    destroy(job_failed);
                end
                % Log task id if failure remains
                if ~isempty(task_failed_remain)
                    [~,FileShort,~] = fileparts(File);
                    [~,SeqFileShort,~] = fileparts(SeqFileName);
                    ErrorFile = fullfile(PathOut,sprintf('%s_%s_errorlog.txt',FileShort,SeqFileShort));
                    fprintf(['Re-run of tasks failed again.\n ... Aborting task and loging failed tasks in ' ErrorFile '\n']);
                    fileID = fopen(ErrorFile, 'at');
                    for tas=1:numel(task_failed_remain)
                        fprintf(fileID,'Task Failed: %d\n', task_failed_remain(tas));
                    end
                    fclose(fileID);
                end
            end
            
            % backup if on
            if Cluster.Backup.On
                % Save temporary for backup
                Backupfile{JobId} = fullfile(Cluster.Backup.Path,sprintf('ResJobs%d.mat',JobId));
                save(Backupfile{JobId},'results','-v7.3');
                % Destroy previous temporary backup
                if JobId>1,delete(Backupfile{JobId-1});end
            end
            % Clean Job
            destroy(job(JobId));
        end
    else
        %% Run local
        for a = 1:NbSimu
            fprintf('Input %d/%d\n',a,NbSimu);
            [Sa, Sphi] = VoxelSim2D_do_one(File,SeqFileName,a,1);
            results{a,1} = Sa;
            results{a,2} = Sphi;
        end
    end
    
    %% Built up the dictionary
    Dico = zeros([size(results,1) size(results{1,1},2)]);
    Error = [];
    for a=1:size(results,1)
        try
            Dico(a,:) = results{a,1} .* exp(1i * results{a,2});
        catch
            Error = [Error a];
        end
    end
    
    %% Save Lookup Table
    if PathOut
        [~,FileShort,~] = fileparts(File);
        [~,SeqFileShort,~] = fileparts(SeqFileName);
        SaveFile = fullfile(PathOut,sprintf('%s_%s_dico.mat',FileShort,SeqFileShort));
        save(SaveFile,'Dico','par','Model','Label','Seq','-v7.3');
    end
end


