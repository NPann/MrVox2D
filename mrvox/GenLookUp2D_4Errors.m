function Dico = GenLookUp2D(ModelFileName,SeqFileName,PathOut,Err)

% This function generate the LookUp table defined in FileIn and save as a file in PathOut.
%
% Input: - FileIn: path to the file where the Model and the MR Sequence are defined
%               (Must be a a cell of strings for multiple files)
%        - PathOut: path to the repertory where the dictionaries will be
%        saved
% Optional: - ServerInfo: Path to the file defining the cluster info for
% parallel jobs (required MCDS + Matlab Runtime Library installed)

% if nargin < 4
    Para = 0;
% elseif nargin == 4
%     Para = 1;
% elseif nargin < 3
%     fprintf('Not enough input argument\nExiting\n');
%     return;
% end

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
%     NbSimu = numel(Err);
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
                if Cluster.MRL.On %Used the standalone version with the MRL routine
                    %                     createTask(job(JobId), @DepVoxelSim2D, 2, {File,Cluster,SimuIdx});
                else % Used the matlab version
                    createTask(job(JobId), @VoxelSim2D_do_one, 2, {File,SeqFileName,SimuIdx,1});
                end
            end
            submit(job(JobId));
            waitForState(job(JobId))
            results(JobSeq(JobId,1):JobSeq(JobId,2),:) = getAllOutputArguments(job(JobId));
            
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
        
        %% Run local
    else
        for a = 1:numel(Err)
            fprintf('Input %d/%d (Err: %d)\n',Err(a),NbSimu,numel(Err));
            [Sa, Sphi] = VoxelSim2D_do_one(File,SeqFileName,Err(a),1);
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
    if ~isempty(Error),fprintf('Error for lines: %d',Error);end
    
    %% Save Lookup Table
%     [~,FileShort,~] = fileparts(File);
%     [~,SeqFileShort,~] = fileparts(SeqFileName);
%     SaveFile = fullfile(PathOut,sprintf('%s_%s_Dico_Err.mat',FileShort,SeqFileShort));
%     save(SaveFile,'Dico','par','Model','Label','Seq','-v7.3');
end


