% Compute feature values for input QSMs for all enabled features.
% Feature can be added with the .APPEND() or .INSERT() methods.
% Feature minimum and maximum values are updated when ever this
% method is called. The returned values are unscaled but can be
% scaled to any interval with the .SCALE() method.
% 
% Inputs:
%
% QSMs          Either a single QSMBCylindrical object or
%               a cell array of such objects.
%
% Outputs:
%
% Value         Matrix of feature values. Each row corresponds to 
%               a single tree, each column to a single feature.
%               The values are not scaled.



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

function Value = compute(ob, QSMs)

    % If only one input QSM, convert to cell array.
    if ~iscell(QSMs)
        QSMs = {QSMs};
    end

    % Number of input trees.
    NTree = numel(QSMs);

    % Check that some trees have been input.
    assert(NTree > 0, 'No models to compute features for.');
    % Check that some features have been set.
    assert(ob.count > 0, 'No features set.');

    % Number of vertical bins.
    NVert = 3;
    % Number of angular bins.
    NAng = 8;
    % Percentage of volume inside vertical bin
    % when determining bin radius.
    VolPer = 0.90;

    % Initialize output matrix of feature values.
    Value = zeros(NTree,ob.count);

    % Go through all trees.
    for iTree = 1:NTree

        % Current tree model.
        QSM = QSMs{iTree};

        % DBH of tree.

        % Try to get diameter from triangulation.
        if QSM.has_property('DBH TRI')
            DBH = QSM.get_property('DBH TRI');

        % Try to get diameter from cylinder fitting.
        elseif QSM.has_property('DBH CYL')
            DBH = QSM.get_property('DBH CYL');

        % By default compute dbh from cylinder model.
        else
            DBH = QSM.dbh();
        end

        % Total model height computed from bounding box.
        TreeHeight = QSM.height();

        % Minimum vertical point of tree.
        TreeMin = QSM.tree_limits(1,3);

        % Total volume of the tree.
        TreeVolume = QSM.volume();

        % Go through all features.
        for iFeature = 1:ob.count

            % Function handle or string identifier of feature.
            fun = ob.handles{iFeature};
            % Feature parameters.
            param = ob.parameters{iFeature};
            % Feature name.
            name = ob.names{iFeature};

            % Function handle features.
            if isa(fun, 'function_handle')

                % Call function with or without parameters.
                if isempty(ob.parameters{iFeature})
                    Val = fun(QSM);
                else
                    Val = fun(QSM, param{:});
                end

            % Built-in features.
            else

                % Each feature is a case in this switch.
                switch fun
                    case 'sb_angle'

                        if not(exist('IStemBranch','var'))

                            % First order branches.
                            IStemBranch = QSM.branch_order == 1;

                        end

                        % Median of stem branch angles.
                        Val = median(QSM.branch_angle(IStemBranch));

                    case 'sb_cluster_size'

                        % Check that enough parameters given.
                        assert( ...
                            numel(param) > 0, ...
                            [ ...
                                'Not enough parameters for feature ''' ...
                                name ...
                                '''.' ...
                            ] ...
                        );

                        % Distance interval from where to lookfor 
                        % cylinders in a cluster.
                        HeightInterval = param{1};

                        if not(exist('IOrd1','var'))

                            % Cylinders that are first in stem 
                            % branches.
                            IOrd1 = find(...
                                QSM.cylinder_branch_order == 1 & ...
                                QSM.cylinder_index_in_branch == 1 ...
                            );
                            %-

                        end

                        % Number of such cylinders.
                        NBranch1 = nnz(IOrd1);

                        % Flag: cylinder has been visited.
                        FVisited = false(NBranch1,1);

                        % Number of found clusters.
                        NClusters = 0;

                        % Sizes of clusters.
                        ClusterSizes = zeros(NBranch1,1);

                        % Average number of 1st order shoots in 
                        % cluster.
                        for iBranch = 1:NBranch1

                            % If visited, i.e., part of some
                            % cluster, skip.
                            if FVisited(iBranch)
                                continue;
                            end

                            % Increase number of clusters.
                            NClusters = NClusters + 1;

                            % Height of new cluster.
                            CurrentHeight = ...
                                QSM.cylinder_start_point(IOrd1(iBranch),3);
                            %-

                            % Logical index of cylinders in
                            % cluster height interval.
                            ILeveled = find(...
                                abs( ...
                                    QSM.cylinder_start_point(IOrd1,3) ...
                                    - CurrentHeight ...
                                ) < HeightInterval & ...
                                not(FVisited)...
                            );

                            % Set found cylinders as visited.
                            FVisited(ILeveled) = true;

                            % Compute cluster size.
                            ClusterSizes(NClusters) = ...
                                length(ILeveled);
                            %-

                        end

                        % Mean value of cluster sizes.
                        Val =  mean(ClusterSizes(1:NClusters));

                    case 'sb_radius'

                        % Check that enough parameters given.
                        assert( ...
                            numel(param) > 0, ...
                            [ ...
                                'Not enough parameters for feature ''' ...
                                name ...
                                '''.' ...
                            ] ...
                        );

                        MaxBranchCount = param{1};

                        if not(exist('IOrd1','var'))

                            % Cylinders that are first in stem
                            % branches.
                            IOrd1 = find(...
                                QSM.cylinder_branch_order == 1 & ...
                                QSM.cylinder_index_in_branch == 1 ...
                            );
                            %-

                        end

                        % Radii of stem branch first cylinders and 
                        % the radii of their stem cylinder parents.
                        Rad1 = [...
                            QSM.cylinder_radius(IOrd1), ...
                            QSM.cylinder_radius( ...
                                QSM.cylinder_parent(IOrd1) ...
                            ) ...
                        ];

                        % Sort from largest to smallest stem
                        % branch radii.
                        [~,ILargest] = sort(Rad1(:,1), 'descend');

                        % Number of ratios to use. At most all of 
                        % the ratios, or the number given as
                        % parameter if it is smaller.
                        NLargest = min( ...
                            MaxBranchCount, ...
                            size(Rad1,1) ...
                        );

                        % Ratio between branch and parent, or if
                        % it is higher than 1 (branches larger 
                        % than parents), set to 1.
                        Val = min(...
                            mean(...
                                Rad1(ILargest(1:NLargest),1) ./ ...
                                Rad1(ILargest(1:NLargest),2)),...
                            1 ...
                        );

                    case 'sb_length'

                        if not(exist('IStemBranch','var'))

                            % First order branches.
                            IStemBranch = QSM.branch_order == 1;

                        end

                        % Ratio between DBH and average stem
                        % branch length.
                        Val = DBH / mean(QSM.branch_length(IStemBranch));


                    case 'sb_distance'

                        % Check that enough parameters given.
                        assert( ...
                            numel(param) > 1, ...
                            [ ...
                                'Not enough parameters for feature ''' ...
                                name ...
                                '''.' ...
                            ] ...
                        );


                        % Width of vertical search window.
                        WindowWidth = param{1};

                        % Maximum value.
                        MaxMeanDist = param{2};

                        if not(exist('IStemBranch','var'))

                            % First order branches.
                            IStemBranch = QSM.branch_order == 1;

                        end

                        % Number of stem branches.
                        NStemBranch = nnz(IStemBranch);

                        % Start heights of stem branches.
                        BranchHeight = ...
                            QSM.branch_height(IStemBranch);
                        %-

                        % Branch distances.
                        VertDists = zeros(NStemBranch^2, 1);

                        % Number of recorded distances.
                        NDist = 0;

                        for iBranch = 1:nnz(IStemBranch)

                            % Logical index of stem branches in 
                            % same vertical bin / window.
                            IBin = BranchHeight < ...
                                BranchHeight(iBranch) + WindowWidth/2 ...
                                & BranchHeight >= ...
                                BranchHeight(iBranch) - WindowWidth/2;
                            %-
                            
                            % Exclude self.
                            IBin(iBranch) = false;

                            % If any branches found in search
                            % window, Compute distances and store.
                            if any(IBin)
                                D = abs(BranchHeight(IBin) ...
                                    - BranchHeight(iBranch));
                                %-

                                % Number of new distances.
                                NBin = size(D,1);

                                % Store new distances.
                                VertDists(NDist+1:NDist+NBin) = D;

                                % Update count.
                                NDist = NDist + NBin;

                            % Otherwise set distance as half of
                            % windown width.
                            else
                                % Store new distance.
                                VertDists(NDist+1) = WindowWidth/2;

                                % Update count.
                                NDist = NDist + 1;
                            end
                        end

                        % Average stem branch distance scaled
                        % by DBH.
                        Val = min(mean(VertDists)/DBH, MaxMeanDist);


                    case 'crown_start_height'

                        if not(exist('ICrown','var'))

                            % Crown cylinders.
                            ICrown = QSM.get_crown();

                        end

                        % Crown cylinders connected to the stem.
                        ICrownStart = ICrown ...
                            & QSM.cylinder_branch_order == 1;
                        %-

                        % Minimum height of starting points.
                        CrownStartHeight = min( ...
                            QSM.cylinder_start_point(ICrownStart,3) ...
                        );

                        % Scale by tree height.
                        Val = (CrownStartHeight ...
                            - QSM.cylinder_start_point(1,3)) ...
                            / TreeHeight;
                        %-

                    case 'crown_height'

                        if not(exist('ICrown','var'))

                            % Crown cylinders.
                            ICrown = QSM.get_crown();

                        end

                        if not(exist('LowestCrownPoint','var'))

                            % Low point of crown cylinders.
                            LowestCrownPoint = min([...
                                QSM.cylinder_start_point(ICrown,3); ...
                                QSM.cylinder_end_point(ICrown,3)...
                            ]);
                        end

                        % Crown total height substracted from
                        % tree height.
                        Val = TreeHeight ...
                            - (LowestCrownPoint ...
                            - QSM.cylinder_start_point(1,3));
                        %-


                    case 'crown_evenness'

                        if not(exist('ICrown','var'))

                            % Crown cylinders.
                            ICrown = QSM.get_crown();

                        end

                        if not(exist('JVertBin','var')) || ...
                            not(exist('JAngBin','var')) || ...
                            not(exist('BinRadii','var'))

                            % Compute vertical and angular bins
                            % and bin radii.
                            [JVertBin, JAngBin, BinRadii] = ...
                                QSMFeatureCalculator.get_cylinder_bins( ...
                                    QSM,...
                                    NVert,... % Num. of vertical bins.
                                    NAng,...  % Num. of angular bins.
                                    VolPer... % Volume percentage.
                            );

                        end

                        % Lower and upper limits of angular bins.
                        CrownAngLimits = [ ...
                            Inf(NAng,1), ...
                            -Inf(NAng,1) ...
                        ];

                        for iAngBin = 1:NAng

                            % Crown cylinders in current bin.
                            IBin = JAngBin == iAngBin & ICrown;

                            if any(IBin)

                                % Minimum height in bin.
                                CrownAngLimits(iAngBin,1) = min([...
                                    QSM.cylinder_start_point(IBin,3); ...
                                    QSM.cylinder_end_point(IBin,3) ...
                                ]);

                                % Maximum height in bin.
                                CrownAngLimits(iAngBin,2) = max([...
                                    QSM.cylinder_start_point(IBin,3); ...
                                    QSM.cylinder_end_point(IBin,3) ...
                                ]);
                            end

                        end

                        % Maximum difference between individual
                        % bin minimum and maximum heights.
                        Val = max(...
                            ( ...
                                min(CrownAngLimits(:,1)) ...
                                - TreeMin ...
                            ) / ( ...
                                max(CrownAngLimits(:,1)) ...
                                - TreeMin ...
                            ),...
                            0 ...
                        );

                    case 'crown_diam_to_height'

                        if not(exist('ICrown','var'))

                            % Crown cylinders.
                            ICrown = QSM.get_crown();

                        end

                        if not(exist('LowestCrownPoint','var'))

                            % Low point of crown cylinders.
                            LowestCrownPoint = min([...
                                QSM.cylinder_start_point(ICrown,3); ...
                                QSM.cylinder_end_point(ICrown,3)...
                            ]);
                        end

                        if not(exist('JVertBin','var')) || ...
                            not(exist('JAngBin','var')) || ...
                            not(exist('BinRadii','var'))

                            % Compute vertical and angular bins
                            % and bin radii.
                            [JVertBin, JAngBin, BinRadii] = ...
                                QSMFeatureCalculator.get_cylinder_bins( ...
                                    QSM,...
                                    NVert,... % Num. of vertical bins.
                                    NAng,...  % Num. of angular bins.
                                    VolPer... % Volume percentage.
                            );

                        end

                        % Height of crown computed as difference
                        % between total tree height and lowest
                        % crown point.
                        CrownTotalHeight = TreeHeight ...
                            - (LowestCrownPoint ...
                            - QSM.cylinder_start_point(1,3));
                        %-

                        % Crown diameter and height ratio.
                        Val = (2*max(BinRadii))/CrownTotalHeight;


                    case 'dbh_to_height'

                        % Check that enough parameters given.
                        assert( ...
                            numel(param) > 0, ...
                            [ ...
                                'Not enough parameters for feature ''' ...
                                name ...
                                '''.' ...
                            ] ...
                        );

                        % Maximum value.
                        MaxVal = param{1};

                        % Ratio between DBH and tree height.
                        % Prevented to be more than given maximum
                        % value.
                        Val = min(DBH / TreeHeight, MaxVal);

                    case 'dbh_to_tree_volume'

                        % Check that enough parameters given.
                        assert( ...
                            numel(param) > 0, ...
                            [ ...
                                'Not enough parameters for feature ''' ...
                                name ...
                                '''.' ...
                            ] ...
                        );

                        % Maximum value.
                        MaxVal = param{1};

                        % Ratio between DBH and tree volume.
                        % Limited from above with the maximum
                        % value.
                        Val = min(DBH/TreeVolume, MaxVal);

                    case 'dbh_to_min_radius'

                        % Check that enough parameters given.
                        assert( ...
                            numel(param) > 0, ...
                            [ ...
                                'Not enough parameters for feature ''' ...
                                name ...
                                '''.' ...
                            ] ...
                        );

                        % Maximum value.
                        MaxVal = param{1};

                        if not(exist('JVertBin','var')) || ...
                            not(exist('JAngBin','var')) || ...
                            not(exist('BinRadii','var'))

                            % Compute vertical and angular bins
                            % and bin radii.
                            [JVertBin, JAngBin, BinRadii] = ...
                                QSMFeatureCalculator.get_cylinder_bins( ...
                                    QSM,...
                                    NVert,... % Num. of vertical bins.
                                    NAng,...  % Num. of angular bins.
                                    VolPer... % Volume percentage.
                            );

                        end

                        % Ratio between DBH and minimum radius
                        % estimate. Capped above with maximum
                        % value.
                        Val = min(DBH / min(BinRadii), MaxVal);

                    case 'volume_below_x_height'

                        % Check that enough parameters given.
                        assert( ...
                            numel(param) > 0, ...
                            [ ...
                                'Not enough parameters for feature ''' ...
                                name ...
                                '''.' ...
                            ] ...
                        );

                        % Height limit below which the volume
                        % percentage is computed.
                        HeightLimit = param{1};

                        % Logical index for stem cylinders.
                        IStem = QSM.cylinder_branch_order == 0;

                        % Vertical volume distribution.
                        % Col 1: cylinder mean point.
                        % Col 2: cylinder volume.
                        VertDist = [...
                            (QSM.cylinder_mid_point(:,3)-TreeMin) ...
                                / TreeHeight,...
                            pi * QSM.cylinder_radius.^2 ...
                                .* QSM.cylinder_length;
                        ];

                        % Branch cylinders below given height.
                        I = VertDist(:,1) < HeightLimit ...
                            & not(IStem);

                        % Sum of selected cylinder volumes divided
                        % by total volume.
                        Val = sum(VertDist(I,2)) / TreeVolume;

                    case 'length_to_volume'

                        % Check that enough parameters given.
                        assert( ...
                            numel(param) > 0, ...
                            [ ...
                                'Not enough parameters for feature ''' ...
                                name ...
                                '''.' ...
                            ] ...
                        );

                        % Maximum value.
                        MaxVal = param{1};

                        % Total branch length divided by total 
                        % tree volume. Capped by max value.
                        Val = min( ...
                            sum(QSM.branch_length)/TreeVolume, ...
                            MaxVal ...
                        );

                    case 'shedding_ratio'

                        if not(exist('JVertBin','var')) || ...
                            not(exist('JAngBin','var')) || ...
                            not(exist('BinRadii','var'))

                            % Compute vertical and angular bins
                            % and bin radii.
                            [JVertBin, JAngBin, BinRadii] = ...
                                QSMFeatureCalculator.get_cylinder_bins( ...
                                    QSM,...
                                    NVert,... % Num. of vert. bins.
                                    NAng,...  % Num. of angu. bins.
                                    VolPer... % Volume percentage.
                            );

                        end

                        % Cylinders is lowest third that are the
                        % first cylinders connected directly to
                        % the stem.
                        ILowBin = JVertBin == 1 ...
                            & QSM.cylinder_branch_order == 1 ...
                            & QSM.cylinder_index_in_branch == 1;
                        %-

                        % Branch indices of such branches.
                        JLowBranch = ...
                            QSM.cylinder_branch_index(ILowBin);
                        %-

                        Val = 1;

                        % If any branches exist in the lowest
                        % third, substract the number of braches
                        % with child branches divided by the total
                        % number of branches.
                        if not(isempty(JLowBranch))
                            Val = Val ...
                                - length(intersect( ...
                                    QSM.branch_parent,JLowBranch ...
                                )) / length(JLowBranch);
                            %-
                        end

                    % Feature id not recognized.
                    otherwise
                        error([ ...
                            'Unable to compute. Built-in feature ' ...
                            'with ID ''' ...
                            fun ...
                            ''' not found.' ...
                        ]);
                end

            end

            % Store individual values into the value matrix.
            Value(iTree, iFeature) = Val;

            % Update current feature min and max values.
            ob.limits(iFeature,1) = min( ...
                Val, ...
                ob.limits(iFeature,1) ...
            );
            ob.limits(iFeature,2) = max( ...
                Val, ...
                ob.limits(iFeature,2) ...
            );

            % Update current feature compute count.
            if isnan(ob.limits(iFeature,3))
                ob.limits(iFeature,3) = 1;
            else
                ob.limits(iFeature,3) = ob.limits(iFeature,3) + 1;
            end

        end
        
        % Clear variables when all features have been computed for
        % a tree.
        clearvars ICrown ...
            JVertBin ...
            JAngBin ...
            BinRadii ...
            IStemBranch ...
            IOrd1 ...
            LowestCrownPoint
        %-

    end

end