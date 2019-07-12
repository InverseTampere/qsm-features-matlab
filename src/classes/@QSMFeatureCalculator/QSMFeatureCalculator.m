% Class for computing QSM feature values. The user can utilize some
% built-in features or append their own and function handle 
% references and optional parameters. Features can be added, listed
% and removed with the respective methods. Default features are 
% available. Once the feature list is final, computations is 
% performed using the .COMPUTE() method. Extreme values of each 
% feature are computed and stored across multiple .COMPUTE() calls,
% and the resulting feature matrix/matrices can be scaled using the
% .SCALE() method.

% This file is part of QSM-Features.
% 
% QSM-Features is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% QSM-Features is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with QSM-Features.  If not, see <http://www.gnu.org/licenses/>.

classdef QSMFeatureCalculator < handle

    properties (Constant)

        % Built-in feature data is set as constant as methods are
        % provided to add and remove built-in features from the 
        % enabled features list. Thus, there should be no need to 
        % manipulate the data directly.

        % String identifiers of features.
        builtin_features = {...
            'sb_angle';...
            'sb_cluster_size';...
            'sb_radius';...
            'sb_length';...
            'sb_distance';...
            'crown_start_height';...
            'crown_height';...
            'crown_evenness';...
            'crown_diam_to_height';...
            'dbh_to_height';...
            'dbh_to_tree_volume';...
            'dbh_to_min_radius';...
            'volume_below_x_height';...
            'length_to_volume';...
            'shedding_ratio';...
        };

        % Feature names.
        builtin_names = {...
            'Stem branch angle';...
            'Stem branch cluster size';...
            'Stem branch radius';...
            'Stem branch length';...
            'Stem branch distance';...
            'Crown start height';...
            'Crown height';...
            'Crown evenness';...
            'Crown diameter / height';...
            'DBH / height ratio';...
            'DBH / tree volume';...
            'DBH / minimum tree radius';...
            'Volume below 55% of height';...
            'Cylinder length / tree volume';...
            'Shedding ratio';...
        };

        % Feature descriptions.
        builtin_desc = {...
            ['Median of the branching angles of the 1st order ' ...
            'branches in degrees. 0 is upwards and 180 downwards. ' ...
            'No parameters.']; ...
            ['Average number of 1st order branches inside a height ' ...
            'interval for 1st order branches. Each branch can only ' ...
            'belong to one interval. Parameter 1: Half of the height ' ...
            'interval size.']; ...
            ['Mean ratio between the N largest 1st order branches ' ...
            'measured at the base and the stem radius at respective ' ...
            'height. Parameter 1: Number of branches (N).']; ...
            ['Average length of 1st order branches normalized by ' ...
            'DBH. No parameters.']; ...
            ['Average distance between 1st order branches computed ' ...
            'using a moving average with a window width. If window ' ...
            'is empty average distance in window is set as half of ' ...
            'window width. Parameter 1: Window width, Parameter 2: ' ...
            'Maximum value for clipping.']; ...
            ['Height of first stem branch in tree crown relative to ' ...
            'tree height. No parameters.']; ...
            ['Vertical distance between the highest and lowest crown ' ...
            'cylinder relative to tree height. No parameters.']; ...
            ['Crown cylinders divided into 8 angular bins. Ratio ' ...
            'between extreme minimum heights in bins. No parameters.']; ...
            'Ratio between crown diameter and height. No parameters.'; ...
            ['Ratio between DBH and total tree height. Parameter 1: ' ...
            'Maximum value.']; ...
            ['Ratio between DBH and total tree volume. Parameter 1: ' ...
            'Maximum value.']; ...
            ['Ratio between DBH and the minimum of the vertical bin ' ...
            'radius estimates. Parameter 1: Maximum value.']; ...
            ['Relative cylinder volume below X% of tree height. ' ...
            'Parameter 1: Relative height limit (X).']; ...
            ['Ratio between total length of all cylinders and total ' ...
            'tree volume. Parameter 1: Maximum value.']; ...
            ['The number of branches without children divided by the ' ...
            'number of all branches in the bottom third. ' ...
            'No parameters.']; ...
        };

        % Default parameters for the built-in functions.
        builtin_parameters = {...
            {};...
            {0.2};...   % Height interval
            {10};...    % Max branch count
            {};...
            {1, 10};... % Window width, max value
            {};...
            {};...
            {};...
            {};...
            {0.05};...  % Max value
            {0.003};... % Max value
            {3};...     % Max value
            {0.55};...  % Height limit
            {10};...    % Max value
            {};...
        };

        % Maximum line width of print-outs.
        line_width = 60;

    end

    properties (Access=private)

        % Function handles for functions computing features.
        handles = {};
        % Cell array of feature names.
        names = {};
        % Cell array of feature descriptions.
        descriptions = {};
        % Cell array of feature parameters.
        % Each element is a cell array of parameters.
        parameters = {};
        % Number of included features.
        count = 0;

        % Maximum length of feature names.
        % Used for printing in ob.list().
        max_name_len = 0;

        % Maximum length of built-in feature ids.
        % Used for printing in ob.list().
        max_id_len = 0;

        % Extreme values of computed features and number of
        % computations. One row per feature. 
        % Columns:
        % 1         Minimum value.
        % 2         Maximum value.
        % 3         Number of computations.
        limits = zeros(0,3);

    end

    methods

        function ob = QSMFeatureCalculator(varargin)
        % Class constructor. 
        %
        % % Initialize with default features.
        % QSMFeatureCalculator()
        % % OR
        % QSMFeatureCalculator('default')
        %
        % % Initialize with no features.
        % QSMFeatureCalculator('none')

            if nargin == 0
                FeatureSet = 'default';
            else
                FeatureSet = varargin{1};
            end

            switch FeatureSet
                case 'default'
                    % By default, add all built-in features.

                    % Number of built-in features.
                    NFeature = numel(ob.builtin_features);

                    % Initialize enabled features properties.
                    ob.handles = cell(NFeature, 1);
                    ob.parameters = cell(NFeature, 1);
                    ob.descriptions = cell(NFeature, 1);
                    ob.names = cell(NFeature, 1);

                    % Copy built-in feature data.
                    for iFeature = 1:NFeature

                        ob.handles{iFeature} = ...
                            ob.builtin_features{iFeature};
                        %-
                        ob.parameters{iFeature} = ...
                            ob.builtin_parameters{iFeature};
                        %-
                        ob.descriptions{iFeature} = ...
                            ob.builtin_desc{iFeature};
                        %-
                        ob.names{iFeature} = ob.builtin_names{iFeature};

                        % Check if new feature name is longer than 
                        % previous maximum.
                        ob.max_name_len = max( ...
                            ob.max_name_len, ...
                            length(ob.names{iFeature}) ...
                        );

                        % Check if new feature identifier is longer than 
                        % previous maximum.
                        ob.max_id_len = max( ...
                            ob.max_id_len, ...
                            length(ob.handles{iFeature}) ...
                        );
                    end

                    % Number of features.
                    ob.count = NFeature;

                case 'none'
                    % Add no features.

                otherwise
                    error(['Unknown feature set: ''' FeatureSet '''']);
                    
            end

            % Initialize feature value upper and lower limits,
            % and computed count.
            ob.limits = nan(ob.count, 3);

        end

        function append(ob, fun, varargin)
        % Append new feature at the end of the feature list.
        % 
        % Examples:
        %
        % FeatCalc.append(Handle, FeatureName);
        % FeatCalc.append(Handle, FeatureName, FeatureDesc);
        % FeatCalc.append(Handle, FeatureName, FeatureDesc, Parameters);
        %
        % Inputs:
        %
        % Handle          Function handle or string identifier of one 
        %                 of the built-in features.
        % FeatureName     Name of the feature. Shown in FeatCalc.list().
        %                 Can be empty if adding a built-in feature.
        % FeatureDesc     String description of the feature. Shown in
        %                 FeatCalc.list('full');
        % Parameters      Optinal parameters to be passed to given function
        %                 handle upon computation or to built-in function.
        %                 Parameters should be a cell array.

            % Position of new feature is at the end of the feature
            % list.
            id = ob.count + 1;

            % Use insert method with new index and input parameters.
            % Insert performs parameter checking.
            ob.insert(id, fun, varargin{:});

        end

        function insert(ob, id, fun, name, desc, param)
        % Insert new feature at position defined by ID number.
        %
        % Examples:
        %
        % FeatCalc.insert(Handle, ID, FeatName);
        % FeatCalc.insert(Handle, ID, FeatName, FeatureDesc);
        % FeatCalc.insert(Handle, ID, FeatName, FeatureDesc, Parameters);
        %
        % Inputs:
        %
        % Handle          Function handle or string identifier of one 
        %                 of the built-in features.
        % ID              Numerical index of where to insert new feature.
        % FeatName        Name of the feature. Shown in FeatCalc.list().
        %                 Can be empty if adding a built-in feature.
        % FeatureDesc     String description of the feature. Shown in
        %                 FeatCalc.list('full');
        % Parameters      Optinal parameters to be passed to given function
        %                 handle upon computation or to built-in function.
        %                 Parameters should be a cell array.

            % Check that index is numeric.
            assert(isnumeric(id),'ID should be numeric.');

            % Check that index within bounds.
            assert( ...
                id > 0 && id <= ob.count + 1, ...
                ['Position ' num2str(id) ' out of bounds.'] ...
            );

            % Default description is empty.
            if nargin < 5
                desc = '';
            end

            % No default parameters.
            if nargin < 6
                param = {};
            end

            % New feature is external.
            if isa(fun,'function_handle')

                % A name should be given if function handle is added.
                assert( ...
                    ischar(name) && ~isempty(name), ...
                    'Second parameter should be a non-empty string.' ...
                );

            % New feature is built-in.
            elseif ischar(fun)

                % Check that string identifier found in built-in options.
                fid = find(ismember(ob.builtin_features, fun));

                if isempty(fid)
                    error([ ...
                        'Built-in feature with ID ''' ...
                        fun ...
                        ''' not found.' ...
                    ]);
                end

                % Use built-in name if not given or empty.
                if nargin < 4 || isempty(name)
                    name = ob.builtin_names{fid};
                end

                % Use built-in description if not given or empty.
                if nargin < 5 || isempty(desc)
                    desc = ob.builtin_desc{fid};
                end

                % Use built-in parameters if not given.
                if nargin < 6
                    param = ob.builtin_parameters{fid};
                else
                    % If given, check that number of parameters matches 
                    % number of default parameters. Check that enough
                    % parameters given.
                    assert( ...
                        numel(param) == ...
                        numel(ob.builtin_parameters{fid}), ...
                        'Parameter count mismatch.' ...
                    );
                end
            else
                error( ...
                    ['First parameter should be a function handle ' ...
                    'or a string identifier of the built-in features.'] ...
                );
            end

            % Update properties appropriately.

            % Add to end (append).
            if id == ob.count + 1

                ob.handles = cat(1, ob.handles, {fun});
                ob.names = cat(1, ob.names, name);
                ob.descriptions = cat(1, ob.descriptions, {desc});
                ob.parameters = cat(1, ob.parameters, {param});
                ob.limits = cat(1, ob.limits, nan(1,3));

            % Add to start.
            elseif id == 1

                ob.handles = cat(1, {fun}, ob.handles);
                ob.names = cat(1, name, ob.names);
                ob.descriptions = cat(1, {desc}, ob.descriptions);
                ob.parameters = cat(1, {param}, ob.parameters);
                ob.limits = cat(1, nan(1,3), ob.limits);

            % Add to middle.
            else

                ob.handles = cat(1,...
                    ob.handles(1:id-1),...
                    {fun},...
                    ob.handles(id:end)...
                );
                ob.names = cat(1,...
                    ob.names(1:id-1),...
                    name,...
                    ob.names(id:end)...
                );
                ob.descriptions = cat(1,...
                    ob.descriptions(1:id-1),...
                    {desc},...
                    ob.descriptions(id:end)...
                );
                ob.parameters = cat(1,...
                    ob.parameters(1:id-1),...
                    {param},...
                    ob.parameters(id:end)...
                );
                ob.limits = cat(1, ...
                    ob.limits(1:id-1,:), ...
                    nan(1,3), ...
                    ob.limits(id:end,:) ...
                );

            end

            % Increase feature count.
            ob.count = ob.count + 1;

        end

        function remove(ob, ids, force)
        % Remove feature with given index.
        %
        % Examples:
        %
        % FeatCalc.remove(ID);
        % FeatCalc.remove(ID, Force);
        %
        % IDs               Numerical indices of feature to be removed.
        % Force             Optinal boolean flag that bypasses 
        %                   confirmation prompt when set to TRUE.

            % Assert legal index value.
            assert(all(ids > 0) && all(ids <= ob.count),...
                   'Feature ID out of range.');
            %-

            % Force flag default if FALSE.
            if nargin < 3
                force = false;
            end

            % Number of ids.
            NID = length(ids);

            for iID = 1:NID

                % Original id.
                ido = ids(iID);

                if NID > 1 && iID > 1

                    id = ido - sum(ids(1:iID-1) < ido);
                else
                    id = ido;
                end

                % Prompt for confirmation if not bypassed.
                if not(force)
                    % Message string.
                    str = sprintf( ...
                        'Removing Feature %d: %s. Continue? (Y/N) [Y]:',...
                        ido,...
                        ob.names{id}...
                    );

                    % Get user input.
                    Button = input(str,'s');

                    % Default value (user presses enter).
                    if isempty(Button)
                        Button = 'y';
                    end

                    % If not 'y', stop function execution.
                    if not(strcmpi(Button,'y'))
                        return;
                    end
                end

                % Update properties appropriately.

                % Remove only feature.
                if ob.count == 1
                    ob.handles = {};
                    ob.names = {};
                    ob.descriptions = {};
                    ob.parameters = {};
                    ob.limits = zeros(0,3);

                % Remove from start.
                elseif id == 1
                    ob.handles = ob.handles(2:end);
                    ob.names = ob.names(2:end);
                    ob.descriptions = ob.descriptions(2:end);
                    ob.parameters = ob.parameters(2:end);
                    ob.limits = ob.limits(2:end,:);

                % Remove last.
                elseif id == ob.count
                    ob.handles = ob.handles(1:end-1);
                    ob.names = ob.names(1:end-1);
                    ob.descriptions = ob.descriptions(1:end-1);
                    ob.parameters = ob.parameters(1:end-1);
                    ob.limits = ob.limits(1:end-1,:);

                % Remove at middle.
                else
                    ob.handles = cat(1,ob.handles(1:id-1), ...
                                       ob.handles(id+1:end));
                    ob.names = cat(1,ob.names(1:id-1), ...
                                     ob.names(id+1:end));
                    ob.descriptions = cat(1,ob.descriptions(1:id-1), ...
                                            ob.descriptions(id+1:end));
                    ob.parameters = cat(1,ob.parameters(1:id-1), ...
                                          ob.parameters(id+1:end));
                    ob.limits = cat(1,ob.limits(1:id-1,:), ...
                                      ob.limits(id+1:end,:));
                    %-

                end
                    
                % Decrease feature count by one.
                ob.count = ob.count - 1;

                % Recompute max length of feature names. Used
                % in .LIST().
                max_name = 0;
                for iFeature = 1:ob.count

                    max_name = max(length(ob.names{iFeature}),max_name);

                end

                ob.max_name_len = max_name;

            end

        end

        function list(ob, varargin)
        % Print a list of feature ID numbers and names of features.
        % Can be used to print enabled features and all available,
        % built-in features. See examples for optional arguments.
        % Argument order and case are ignored.
        %
        % Examples:
        %
        % % Print enabled features in short form.
        % FeatCalc.list();
        %
        % % Print enabled features in long form, i.e.,
        % % include feature descriptions.
        % FeatCalc.list('full');
        % % OR
        % FeatCalc.list('verbose');
        %
        % % Print available features in short form.
        % FeatCalc.list('builtin');
        % % OR
        % FeatCalc.list('built-in');
        % % OR
        % FeatCalc.list('available');
        %
        % % Print available features in long form, i.e.,
        % % include feature descriptions.
        % FeatCalc.list('builtin','verbose');
        % % OR
        % FeatCalc.list('built-in','full');
        % % OR
        % FeatCalc.list('available','verbose');
        %

            % Flag for including description when TRUE.
            FDesc = false;
            % Flag for using built-in features instead of 
            % enabled features. List built-in when TRUE.
            FBuiltin = false;

            % Parse optional arguments.
            if nargin > 1

                for iArgin = 1:numel(varargin)

                    par = varargin{iArgin};

                    switch lower(par)

                        case {'full','verbose'}
                            FDesc = true;
                        case {'builtin','built-in','available'}
                            FBuiltin = true;
                        otherwise
                            warning([ ...
                                'Ignoring unknown input parameter: ''' ...
                                par ...
                                '''' ...
                            ]);
                    end

                end

            end

            % List built-in features.
            if FBuiltin
                % Get feature properties.
                IDs = ob.builtin_features;
                Names = ob.builtin_names;
                Descs = ob.builtin_desc;

                % Print header with count.
                fprintf( ...
                    '\nAvailable built-in features (%d)\n', ...
                    numel(ob.builtin_features) ...
                );

                % Print feature list.
                ob.print_features(IDs, Names, Descs, FDesc);

            % List enabled features.
            else
                % Check that any features are present.
                assert(ob.count > 0, 'No features set.');

                % Get feature properties.
                IDs = num2cell(1:ob.count);
                Names = ob.names;
                Descs = ob.descriptions;

                % Print header with count.
                fprintf('\nEnabled features (%d)\n', ob.count);
                % Print feature list.
                ob.print_features(IDs, Names, Descs, FDesc);
            end

        end

        function init_limits(ob)
        % Initialize feature minumum and maximum values.

            % Set limits as NaNs.
            ob.limits = nan(size(ob.limits));

        end

        function Lim = get_limits(ob)
        % Return feature minimum and maximum values accumulated
        % during .COMPUTE() calls.

            % Return limits.
            Lim = ob.limits(:,1:2);

        end

        function Scaled = scale(ob, Val, Lim)
        % Scale feature values to given limits based on stored
        % feature minimum and maximum values. These values are updated
        % when .COMPUTE() is used. If the minimum and maximum values of
        % a given feature are equal, all the values are mapped onto the
        % lower limit.
        % 
        % Inputs:
        %
        % Val           Feature value matrix. Each row corresponds to 
        %               a single tree, each column to a single feature.
        %               Column count should match feature count.
        % Lim           Range to which the values are mapped to. By 
        %               default the range is [0, 1]. Value should be a
        %               two-element vector [min, max].
        %
        % Examples:
        %
        % Scaled = FeatCalc.scale(Val);
        % Scaled = FeatCalc.scale(Val, Lim);
        % Scaled = FeatCalc.scale(Val, [-1, 1]);

            % Check that column count matches feature count.
            assert( ...
                size(Val,2) == size(ob.limits,1), ...
                ['Feature count mismatch. Feature count probably ' ...
                'changed after running COMPUTE.'] ...
            );

            % Check that only numeric values present in limits.
            assert( ...
                not(any(isnan(ob.limits(:)))), ...
                'NaNs present in limits.' ...
            );

            % Default limits.
            if nargin < 3
                Lim = [0, 1];
            end

            % Check that limit vector has correct length.
            assert( ...
                length(Lim) == 2, ...
                'Limit vector should have length of 2.' ...
            );

            % Check that lower limit smaller than upper limit.
            assert( ...
                Lim(2) > Lim(1), ...
                'Lower limit above upper limit.' ...
            );

            % Initialize scaled values.
            Scaled = zeros(size(Val));

            % Number of features.
            NFeature = size(Val,2);
            % Number of trees.
            NTree = size(Val,1);

            % Go through all features.
            for iFeature = 1:NFeature
                % Get minimum and maximum values stored in the object.
                mi = ob.limits(iFeature,1);
                ma = ob.limits(iFeature,2);

                % If min and max are the same, avoid NaN by setting
                % divider as one rather than zero.
                if ma == mi
                    d = 1;
                else
                    d = ma - mi;
                end

                % Go through all the trees.
                for iTree = 1:NTree

                    % Short format for feature value.
                    v = Val(iTree, iFeature);

                    % Linear interpolation between stored min and max 
                    % values. Furher scaled to given limits.
                    Scaled(iTree, iFeature) = (v - mi)/d ...
                        *(Lim(2)-Lim(1)) + Lim(1);
                    %-

                end

            end

        end

        % External definition.
        Value = compute(ob, QSMs)

    end

    % Internal helper functions.
    methods (Access=private)

        function print_features(ob, IDs, Names, Descs, FDesc)
        % Function to handle feature list printing.
        %
        % Inputs:
        %
        % IDs           Cell-array of numeric or string identifiers.
        % Names         Cell-array of feature name strings.
        % Descs         Cell-array of feature descriptions strings.
        % FDesc         Flag which includes descriptions when TRUE.

            % Number of features to print.
            NFeature = numel(Names);

            % String to print between columns.
            ColBreak = '   ';

            % Check if IDs are numeric or strings and set 
            % print format and column with appropriately.
            if isnumeric(IDs{end})
                fmt = 'd';
                NDigit = max(2,length(num2str(ob.count)));
            else
                fmt = 's';
                NDigit = max(2,ob.max_id_len);
            end

            % Number of character to print on limiting lines.
            NDash = NDigit+length(ColBreak)+ob.max_name_len;

            % Maximum width of printed line.
            % Constant maximum or le
            MaxWidth = max(ob.line_width, NDash);

            % If descriptions are included, use maximum width 
            % for limiting lines.
            if FDesc
                NDash = MaxWidth;
            end

            fprintf('\n');
            % Top rule.
            fprintf(repmat('=',1,NDash));
            fprintf('\n');
            % Column headers.
            fprintf(['%-' num2str(NDigit) 's'], 'ID');
            fprintf(ColBreak);
            fprintf('%s', 'Name');
            fprintf('\n');

            % Mid rule.             
            fprintf(repmat('–',1,NDash));
            fprintf('\n');

            % Go through all features.
            for iFeature = 1:NFeature

                % Print ID.
                fprintf(['%-' num2str(NDigit) fmt], IDs{iFeature});
                fprintf(ColBreak);
                % Print fearure name.
                fprintf('%s', Names{iFeature});
                fprintf('\n');

                % If description is included, print description
                % word-by-word so maximum width is not exceeded.
                if FDesc

                    % Size of indent before description.
                    NIndent = NDigit + length(ColBreak);

                    % Get individual words separated by spaces.
                    DescWords = strsplit(Descs{iFeature}, ' ');

                    % Number of words in description.
                    NWord = numel(DescWords);

                    % Width of available line after indent.
                    LineWidth = MaxWidth - NIndent;

                    % Number of characters left of current line.
                    LineLeft = LineWidth;

                    % Store indent that is printed before each description
                    % line.
                    Indent = sprintf( ...
                        ['%-' num2str(NDigit) 's' ColBreak], ...
                        '|' ...
                    );

                    % Print first line indent.
                    fprintf('%s', Indent);

                    % Go through all the words.
                    for iWord = 1:NWord

                        % Length of next word.
                        NextLength = length(DescWords{iWord});

                        % Compute how many characters would
                        % be left if word printed.
                        LineLeft = LineLeft - NextLength;

                        % If word does no fit, start new line.
                        if LineLeft < 0

                            % Print line indent.
                            fprintf('\n%s', Indent);
                            % Initialize line width.
                            LineLeft = LineWidth - NextLength;
                        end

                        % Print word.
                        fprintf('%s', DescWords{iWord});

                        % If the word is not the last one, print 
                        % trailing space, if it fits the line.
                        % Otherwise next word starts a new line
                        % and no space is required.
                        if LineLeft > 0 && iWord < NWord
                            fprintf(' ');
                            % Update available line length.
                            LineLeft = LineLeft - 1;
                        end

                    end

                    % Empty line after description.
                    fprintf('\n\n');
                end

            end
            % Bottom rule.
            fprintf(repmat('=',1,NDash));
            fprintf('\n');

        end

    end

    methods (Static)

        function [JVertBin, JAngBin, BinRadii] = get_cylinder_bins( ...
            QSM, ...
            NVert, ...
            NAng, ...
            VolPer ...
        )
        % Divide the QSM cylinders into vertical and angular bins.
        %
        % Inputs:
        %
        % QSM           QSMBCylindrical object that contains the cylinders.
        % NVert         Number of vertical bins.
        % NAng          Number of angular bins.
        % VolPer        Percentage of volume required to define vertical
        %               bin radius.
        %
        % Outputs:
        %
        % JVertBin      Index of vertical bin of each cylinder.
        % JAngBin       Index of angular bin or each cylinder.
        % BinRadii      Vector of radii of vertical bins.

            % Minimum vertical point of tree.
            TreeMin = QSM.tree_limits(1,3);

            % Total tree height.
            TreeHeight = QSM.height();

            % Volume of each cylinder.
            CylinderVolume = pi * QSM.cylinder_radius.^2 ...
                .* QSM.cylinder_length;

            % Cylinder mid points.
            CylinderMean = QSM.cylinder_mid_point;

            % Logical index for stem cylinders.
            IStem = QSM.cylinder_branch_order == 0;
            % And other cylinders.
            INotStem = not(IStem);

            % Number of cylinders.
            NCyl = QSM.block_count;

            % Relative vertical bin limits.
            VertBinLimits = (1:NVert)'/NVert;
            VertBinLimits = [[0; VertBinLimits(1:end-1)],VertBinLimits];

            % Angular bin limits.
            AngBinLimits = linspace(-pi,pi,NAng+1);
            AngBinLimits = [AngBinLimits(1:end-1); AngBinLimits(2:end)];

            % Index of the vertical bin of each cylinder. 1-indexed.
            JAngBin = zeros(NCyl,1);

            % Index of the angular bin of each cylinder. 1-indexed.
            JVertBin = zeros(NCyl, 1);

            % Vertical distribution based on cylinder center points.
            VertDist = (CylinderMean(:,3) - TreeMin) / TreeHeight;

            BinRadii = zeros(NVert, 1);
            BinCenters = zeros(NVert, 2);

            % Iterate over vertical bins.
            % - Find indices of cylinders in each vertical bin.
            % - Find horizontal mean point in each vertical bin and
            %   assign bin cylinders into angular bins around that 
            %   mean point.
            for iBin = 1:NVert

                % Indices of cylinders in current vertical bin.
                switch iBin
                    case 1
                        IBin = VertDist(:,1) <  VertBinLimits(iBin,2);

                    case NVert
                        IBin = VertDist(:,1) >= VertBinLimits(iBin,1);

                    otherwise
                        IBin = VertDist(:,1) >= VertBinLimits(iBin,1) & ...
                               VertDist(:,1) <  VertBinLimits(iBin,2);
                    %-

                end

                % Store bin information.
                JVertBin(IBin) = iBin;

                % Bin stem cylinders.
                IBinStem   = IBin & IStem;
                % Bin branch cylinders.
                IBinBranch = IBin & INotStem;

                % If bin branches exist.
                if nnz(IBinBranch)

                    % Compute an estimate of bin radius.

                    % Volume of branches in bin.
                    BinBranchVol = sum(CylinderVolume(IBinBranch));

                    % If bin contains stem, compute mean in XY-plane.
                    if nnz(IBinStem)
                        BinCenters(iBin,:) = ...
                            mean(CylinderMean(IBinStem,1:2),1);
                        %-

                    else
                        % No stem present in vertical bin, use stem mean
                        % point from lower vertical bin, or use branch
                        % cylinder mean if no lower bin exist.
                        if iBin > 1
                            BinCenters(iBin,:) = BinCenters(iBin-1,:);
                        else
                            BinCenters(iBin,:) = ...
                                mean(CylinderMean(IBinBranch,1:2),1);
                            %-
                        end
                    end

                    % Radial volume distribution of bin: cylinder volume
                    % as a  function of horizontal distance from stem
                    % mean. Sorted by distance.
                    RVD = [...
                        sqrt(sum(bsxfun(...
                            @minus,...
                            CylinderMean(IBinBranch,1:2),...
                            BinCenters(iBin,:)).^2,2)...
                        ),...
                        CylinderVolume(IBinBranch)
                    ];

                    % Sort from closest to farthest from center.
                    [~,IBinRD] = sort(RVD(:,1),'ascend');

                    % Find distance at which <VolPer> percentage of
                    % volume is reached. Store index of such cylinder.
                    JRadLim = find( ...
                        cumsum(RVD(IBinRD,2)/BinBranchVol) > VolPer, ...
                        1, ...
                        'first' ...
                    );

                    % Store distance as radius estimate.
                    BinRadii(iBin) = RVD(IBinRD(JRadLim),1);

                % Vertical bin does not contain branch cylinders.
                else

                    % Set radii as mean radius of stem cylinders in bin.
                    BinRadii(iBin) = mean(QSM.cylinder_radius(IBinStem));
                end

                % Angular distribution.

                % Location in 2D of branch cylinders relative to bin
                % center.
                BinCylRelLoc = bsxfun(...
                    @minus,...
                    CylinderMean(IBinBranch,1:2),...
                    BinCenters(iBin,:)...
                );

                % Transform locations into polar coordinates and store
                % angle.
                [BinTheta, ~] = cart2pol( ...
                    BinCylRelLoc(:,1), ...
                    BinCylRelLoc(:,2) ...
                );

                % Convert logical index into numerical.
                JBinBranch = find(IBinBranch);

                for iAngBin = 1:NAng

                    % Find cylinders in each bin based on angle.
                    IAngBin = BinTheta > AngBinLimits(1,iAngBin) ...
                            & BinTheta <= AngBinLimits(2,iAngBin);

                    % Store respective angular bin index.
                    JAngBin(JBinBranch(IAngBin)) = iAngBin;
                end
            end

        end

    end

end