

module load matlab/r2022b

matlab -singleCompThread -batch "addpath(genpath('~/Documents/MATLAB/')); rng('shuffle'); DA_regression_submit_batch;  exit"