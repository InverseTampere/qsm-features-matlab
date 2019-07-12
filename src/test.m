%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Add QSM-Blocks classes from where you have them.  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add class definitions to path.
addpath('./classes');

% Initialize the calculator.
FC = QSMFeatureCalculator();

% List currently activate features. 
% Built-in features are activate by default.
FC.list();

% Remove features 1, 4 and 15. Confirmation is asked for each.
% Confirmation can be bypassed by setting the optional second
% argument to TRUE, i.e., FC.remove([1,4,15],true).
FC.remove([1,4,15]);

% Check changes.
FC.list();

% List built-in features that can be activated.
FC.list('builtin');

% Add Shedding ratio feature back in, i.e., activate it.
% The string identifier for the built-in features are listed 
% in the .list('builtin') output.
FC.append('shedding_ratio');

% Create an inline function to add a custom feature.
% Any function handle would work the same way.
% Function should have a QSM as the first argument
% (optional arguments possible), and return a scaler
% feature value.
x = @(QSM) mean(QSM.cylinder_radius./QSM.cylinder_length);

% Add new feature defined by inline function.
% This time INSERT is used instead of APPEND.
% The syntax is the same, other than a numeric
% index is also given. Here we append the new
% feature as Feature 1.
%
% Feature name and description are also given.
% Optional arguments can also be passed to the 
% custom feature calculation function. See 
% QSMFeatureCalculator.append for details.
FC.insert( ...
    1, ...
    x, ...
    'Rad-Len ratio', ...
    'Average ratio between cylinder length and radius' ...
);

% List activated features with their descriptions.
% Shedding ratio has been added to the bottom of the 
% table, and Rad-Len ratio to the top.
FC.list('full');

% Create simple example QSM.
QSM = QSMBCylindrical('example');

% Compute feature values for the single QSM.
FeatVal = FC.compute(QSM);

% Show results.
disp('Computed feature values for a single QSM.');
disp(FeatVal(:));

% Values can be scaled and mapped on to a specific range.
% However as the values were computed for a single QSM,
% the extreme values of each feature are equal. Thus all
% the values will be mapped to the lower limit (-1).
Scaled = FC.scale(FeatVal,[-1 1]);

% Show scaled values.
disp('Scaled feature values.');
disp(Scaled(:));

%%
% Check recorded extreme values. One row per feature,
% min value in first column, max in second.
Lim = FC.get_limits()

% Limits can be initialized with the init_limits method.
FC.init_limits();

% Re-check limits. All have been reset to NaNs.
FC.get_limits()