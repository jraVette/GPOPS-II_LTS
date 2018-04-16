%The purpose of this DOE is to get parallel jobs running in true Parallel
%in Palmetto.
%Updated 24 Jul 2017 - Jeff Anderson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c


%Flags and constants
nRepeats = 10;                                                              %How many of each iterate to run
setupEachIterateDirectory   = true;                                        %true or false, true will create a folder in the ./currentlyRunning folder for each iterate, false will just setup the DAQ files and create a summary results .mat file
currentlyRunningDirectory   = 'currentlyRunning';                          %Name of the currently running directory

daqFilenameToSave           = 'daqFile.mat';                               %what name to save the daqFile (should probably use daqFile.mat)
resultsFilename             = 'results.mat';
compiledVersion             = false;
templateIterateDirectory    = 'templateIterate';                           %What folder is the template iterate
matlabCodeDirectory         = 'compilingSimulation';




%The DOE cell array is where to make the applicable runs. This cell array
%is a row vector. Each row is a new set of experiments with columns fields
%to change and values to chage to. note that values column can have
%multiple columns but must be consistent as it will use each column to
%setup the experiment.

%REV01 Course DOE                  1        2       3      4  5    6    7      8
doe = {...
    {'header.horizon'             [200 300 350 400 450 500 500 500 500 500 550 600 700 800 ] 
     'header.controlHorizon'      [50  50  50  50  50  30  40  50  60  100 50  50  50  50]}; %Horizon DOE
};


%REV03 Repeat 
doe = {...
    {'header.horizon'             [300] 
     'header.controlHorizon'      [60]}; %Horizon DOE
};

%REV04 Repeat w/ lower mesh refinement
doe = {...
    {'header.horizon'             [300] 
     'header.controlHorizon'      [60]
     'header.setup.mesh.maxiterations', [1]};
};

%REV05 Col points DOE
doe = {...
    {'header.horizon'             		[300   300   300] 
     'header.controlHorizon'      		[60    60    60 ]
     'header.setup.mesh.maxiterations', [2     4     2  ]
     'header.setup.colpointsmin'        [2     2     2  ]
     'header.setup.colpointsmax'        [6     6     2  ]};
};

%REV06 More col point and look at scaling
doe = {...
    {'header.setup.mesh.maxiterations', [2     4     6   8 10  ]};
    {'header.setup.colpointsmin'        [1     1     1   1      2     2  2 ]
     'header.setup.colpointsmax'        [1     2     10   4      2     3  10 ]};
};

doe = {...
    {'header.setup.scales.method'       { 'automatic-bounds', 'automatic-guess', 'automatic-guessUpdate', 'automatic-hybridUpdate'}};
};

%REV07 I setup the col points wrong last time, so look at them
doe = {...
     {'header.setup.mesh.colpointsmin'        [1     1     1   1      2     2  2 ]
      'header.setup.mesh.colpointsmax'        [1     2     10   4      2     3  10 ]};
};

%REV08 vx Let's see
doe = {...
    {'header.horizon'             [300] 
     'header.controlHorizon'      [60]}; %Horizon DOE
}

        
%Setup each daq DOE
addpath('compilingSimulation/')
fprintf('Number of master experiments = %i\n',numel(doe))
count = 1;
for i = 1:numel(doe)
    
    localDoe = doe{i};
    values = localDoe{1,2};
    daq = generateInitialDaq;
    for j = 1:numel(values)
        
        for iField = 1:size(localDoe,1)
            key = localDoe{iField,1};
            value = localDoe{iField,2}(j);  
            if isa(value,'cell')
                temp = value;
                clear value
                value = temp{1};
                value = ['''' value ''''];
            else
                value = num2str(value);
            end
            evalc(['daq.' key '=' value ]);
            daq.header.doe.keys{iField} = localDoe{iField,1};
            if isa(localDoe{iField,2},'cell')
                daq.header.doe.values{iField} = localDoe{iField,2}{j};
            else
                daq.header.doe.values(iField) = localDoe{iField,2}(j);
            end
        end       
        
        filename{count} = sprintf('DOE%03i_%03i',i,j);
        daq.header.doe.name = filename{count};
        daq.header.filename = daqFilenameToSave;
        DOE_DAQ{count} = daq;
        count = count+1;
    end
end
fprintf('Number of local experiments = %i\n',count)


%Copy files from template iterate and save daq files in correct directory over so that we can run simulation
if setupEachIterateDirectory
    for iFile = 1:length(filename)
        for iRepeat = 1:nRepeats
            if nRepeats > 1
                currentFile = fullfile(currentlyRunningDirectory,[filename{iFile} sprintf('_Repeat%03i',iRepeat)]);
            else
                currentFile = fullfile(currentlyRunningDirectory,filename{iFile});
            end
            mkdir(currentFile)
            if compiledVersion
                sysCommand = sprintf('cp -r %s/* %s',templateIterateDirectory,currentFile);
                system(sysCommand);
            else
                sysCommand = sprintf('cp %s/runBatchMpcSimulationMatlabVersion.m %s',matlabCodeDirectory,currentFile);
                system(sysCommand)
            end
            
            currentDirectory = pwd;
            cd(currentFile);

            daq = DOE_DAQ{iFile};
            daq.header.path = pwd;
            save(daqFilenameToSave,'daq')
            cd(currentDirectory);
        end
    end
end



    
    
    
