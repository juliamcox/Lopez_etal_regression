

% Get a handle to the cluster
c = parcluster;

%% Required arguments in order to submit MATLAB job

% Specify the walltime (e.g. 4 hours)
c.AdditionalProperties.WallTime = '48:00:00';

% Specify an account to use for MATLAB jobs (e.g. pXXXX, bXXXX, etc)
c.AdditionalProperties.AccountName = 'p31438';

% Specify a queue/partition to use for MATLAB jobs (e.g. short, normal, long)
c.AdditionalProperties.QueueName = 'normal';

% Constrain this MPI job
%c.AdditionalProperties.Constraint = '"[quest9|quest10|quest11|quest12]"';

%% optional arguments but worth considering

% Specify memory to use for MATLAB jobs, per core (default: 4gb)
c.AdditionalProperties.MemUsage = '10gb';

% Specify number of nodes to use
c.AdditionalProperties.Nodes = 1;

% Specify e-mail address to receive notifications about your job
c.AdditionalProperties.EmailAddress = 'julia.cox@northwestern.edu';

% The script that you want to run through SLURM needs to be in the MATLAB PATH
% Here we assume that quest_parallel_example.m lives in the same folder as submit_matlab_job.m
addpath(pwd)

% Finally we will submit the MATLAB script quest_parallel_example to SLURM such that MATLAB
% will request enough resources to run a parallel pool of size 4 (i.e. parallelize over 4 CPUs).,
job = c.batch('DA_activeAvoid_script', 'Pool', 16, 'CurrentFolder', '.');