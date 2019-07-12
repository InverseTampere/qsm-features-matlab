# QSM-Features

Quantitative structure models - Features is a Matlab class definition for computing tree model features for applications such as tree species classification. Features can utilize both geometric and topological information from QSMs.

## Description

The repository provides the `QSMFeatureCalculator` class that contains holds a list of feature definitions. Value of each feature can be computed for either a single QSM or a set of QSMs. Features described in the publication [Åkerblom et al.: Non-intersecting leaf insertion algorithm for tree structure models](http://rsfs.royalsocietypublishing.org/content/8/2/20170045) are built-in and their implementations are part of the source. Additionally, users can add their own features by providing the handle of a custom function that receives a single QSM and returns a scalar feature value. Features can be listed, added, removed. 

The `compute`-method receives the QSMs for which values should be computed for and returns the computed feature values for each QSM and each feature as a matrix, where each row corresponds to a single QSM and each column to a single feature. 

The `QSMFeatureCalculator` object stores the extreme values of each feature across multiple `compute`-calls. This so that the feature values can be later scaled with the extreme values, using the `scale`-method.

## Dependencies

### QSMB and QSMBCylindrical

Quantitative structure models (QSMs) are expected to be defined as objects of QSMB subclasses, such as QSMBCylindrical. These Matlab classes are part of the [QSM-Blocks](https://github.com/InverseTampere/qsm-blocks-matlab) repository.

## Basic usage

### Initialization

Initialize the calculator.

```Matlab
FC = QSMFeatureCalculator();
```

### Listing, adding and removing features

List currently activate features. 
Built-in features are activate by default.

```Matlab
FC.list();

% Enabled features (15)
% 
% ==================================
% ID   Name
% ––––––––––––––––––––––––––––––––––
% 1    Stem branch angle
% 2    Stem branch cluster size
% 3    Stem branch radius
% 4    Stem branch length
% 5    Stem branch distance
% 6    Crown start height
% 7    Crown height
% 8    Crown evenness
% 9    Crown diameter / height
% 10   DBH / height ratio
% 11   DBH / tree volume
% 12   DBH / minimum tree radius
% 13   Volume below 55% of height
% 14   Cylinder length / tree volume
% 15   Shedding ratio
% ==================================
```

Remove features 1, 4 and 15. Confirmation is asked for each.
Confirmation can be bypassed by setting the optional second
argument to TRUE, i.e., `FC.remove([1,4,15],true)`.

```Matlab
FC.remove([1,4,15]);
```

Check changes.

```Matlab
FC.list();

% Enabled features (12)
% 
% ==================================
% ID   Name
% ––––––––––––––––––––––––––––––––––
% 1    Stem branch cluster size
% 2    Stem branch radius
% 3    Stem branch distance
% 4    Crown start height
% 5    Crown height
% 6    Crown evenness
% 7    Crown diameter / height
% 8    DBH / height ratio
% 9    DBH / tree volume
% 10   DBH / minimum tree radius
% 11   Volume below 55% of height
% 12   Cylinder length / tree volume
% ==================================
```

List built-in features that can be activated.

```Matlab
FC.list('builtin');

% Available built-in features (15)
% 
% =====================================================
% ID                      Name
% –––––––––––––––––––––––––––––––––––––––––––––––––––––
% sb_angle                Stem branch angle
% sb_cluster_size         Stem branch cluster size
% sb_radius               Stem branch radius
% sb_length               Stem branch length
% sb_distance             Stem branch distance
% crown_start_height      Crown start height
% crown_height            Crown height
% crown_evenness          Crown evenness
% crown_diam_to_height    Crown diameter / height
% dbh_to_height           DBH / height ratio
% dbh_to_tree_volume      DBH / tree volume
% dbh_to_min_radius       DBH / minimum tree radius
% volume_below_x_height   Volume below 55% of height
% length_to_volume        Cylinder length / tree volume
% shedding_ratio          Shedding ratio
% =====================================================
```

Add Shedding ratio feature back in, i.e., activate it.
The string identifier for the built-in features are listed 
in the .list('builtin') output.

```Matlab
FC.append('shedding_ratio');
```

Create an inline function to add a custom feature.
Any function handle would work the same way.
Function should have a QSM as the first argument
(optional arguments possible), and return a scaler
feature value.

```Matlab
x = @(QSM) mean(QSM.cylinder_radius./QSM.cylinder_length);
```

Add new feature defined by inline function.
This time `insert` is used instead of `append`.
The syntax is the same, other than a numeric
index is also given. Here we append the new
feature as Feature 1.

Feature name and description are also given.
Optional arguments can also be passed to the 
custom feature calculation function. See 
`QSMFeatureCalculator.append` for details.

```Matlab
FC.insert( ...
    1, ...
    x, ...
    'Rad-Len ratio', ...
    'Average ratio between cylinder length and radius' ...
);
```

List activated features with their descriptions.
*Shedding ratio* has been added to the bottom of the 
table, and *Rad-Len* ratio to the top.

```Matlab
FC.list('full');

% Enabled features (14)
% 
% ============================================================
% ID   Name
% ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
% 1    Rad-Len ratio
% |    Average ratio between cylinder length and radius
% 
% 2    Stem branch cluster size
% |    Average number of 1st order branches inside a height 
% |    interval for 1st order branches. Each branch can only 
% |    belong to one interval. Parameter 1: Half of the height
% |    interval size.
% 
% 3    Stem branch radius
% |    Mean ratio between the N largest 1st order branches 
% |    measured at the base and the stem radius at respective 
% |    height. Parameter 1: Number of branches (N).
% 
% 4    Stem branch distance
% |    Average distance between 1st order branches computed 
% |    using a moving average with a window width. If window 
% |    is empty average distance in window is set as half of 
% |    window width. Parameter 1: Window width, Parameter 2: 
% |    Maximum value for clipping.
% 
% 5    Crown start height
% |    Height of first stem branch in tree crown relative to 
% |    tree height. No parameters.
% 
% 6    Crown height
% |    Vertical distance between the highest and lowest crown 
% |    cylinder relative to tree height. No parameters.
% 
% 7    Crown evenness
% |    Crown cylinders divided into 8 angular bins. Ratio 
% |    between extreme minimum heights in bins. No parameters.
% 
% 8    Crown diameter / height
% |    Ratio between crown diameter and height. No parameters.
% 
% 9    DBH / height ratio
% |    Ratio between DBH and total tree height. Parameter 1: 
% |    Maximum value.
% 
% 10   DBH / tree volume
% |    Ratio between DBH and total tree volume. Parameter 1: 
% |    Maximum value.
% 
% 11   DBH / minimum tree radius
% |    Ratio between DBH and the minimum of the vertical bin 
% |    radius estimates. Parameter 1: Maximum value.
% 
% 12   Volume below 55% of height
% |    Relative cylinder volume below X% of tree height. 
% |    Parameter 1: Relative height limit (X).
% 
% 13   Cylinder length / tree volume
% |    Ratio between total length of all cylinders and total 
% |    tree volume. Parameter 1: Maximum value.
% 
% 14   Shedding ratio
% |    The number of branches without children divided by the 
% |    number of all branches in the bottom third. No 
% |    parameters.
% 
% ============================================================
```

### Computing values

Create simple example QSM.

```Matlab
QSM = QSMBCylindrical('example');
```

Compute feature values for the single QSM.

```Matlab
FeatVal = FC.compute(QSM);
```

Show results.

```Matlab
disp('Computed feature values for a single QSM.');
disp(FeatVal(:));

% Computed feature values for a single QSM.
%     0.1406
%     3.0000
%     0.6667
%          0
%     0.5000
%     1.0000
%          0
%     6.0000
%     0.0500
%     0.0030
%     1.3333
%          0
%     9.5214
%     1.0000
```

### Scaling values

Values can be scaled and mapped on to a specific range.
However as the values were computed for a single QSM,
the extreme values of each feature are equal. Thus all
the values will be mapped to the lower limit (-1).

```Matlab
Scaled = FC.scale(FeatVal,[-1 1]);
```

Show scaled values.

```Matlab
disp('Scaled feature values.');
disp(Scaled(:));

% Scaled feature values.
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
%     -1
```

Check recorded extreme values. One row per feature,
min value in first column, max in second.

```Matlab
Lim = FC.get_limits()

% Lim =
% 
%     0.1406    0.1406
%     3.0000    3.0000
%     0.6667    0.6667
%          0         0
%     0.5000    0.5000
%     1.0000    1.0000
%          0         0
%     6.0000    6.0000
%     0.0500    0.0500
%     0.0030    0.0030
%     1.3333    1.3333
%          0         0
%     9.5214    9.5214
%     1.0000    1.0000
```

Limits can be initialized with the `init_limits` method.

```Matlab
FC.init_limits();
```

Re-check limits. All have been reset to NaNs.

```Matlab
FC.get_limits()

% ans =
%
%    NaN   NaN
%    NaN   NaN
%    NaN   NaN
%    NaN   NaN
%    NaN   NaN
%    NaN   NaN
%    NaN   NaN
%    NaN   NaN
%    NaN   NaN
%    NaN   NaN
%    NaN   NaN
%    NaN   NaN
%    NaN   NaN
%    NaN   NaN
```
