function [out_featuresSyn, out_labelsSyn] = SMOTE(in_features, in_labels, in_k, in_SMOTE)

% this function implements the SMOTE  method as proposed in the following paper:
%
% N. Chawla, K. Bowyer, L. Hall, and W. Kegelmeyer. Smote: synthetic minority 
% over-sampling technique. Arxiv preprint arXiv:1106.1813, 2011.
%
% the purpose of the ADASYN method is to improve class balance towards
% equally-sized classes for a given input dataset. this is achieved by
% synthetically creating new examples from the minority class via linear
% interpolation between existing minority class samples. this approach is
% known as the SMOTE method
%
%
% INPUTS:
%----------
%in_features:
%(N \times P) matrix of numerical features. each row is one example, each
%column is one feature, hence there are N examples with P features each.
%
%in_labels:
%boolean N-vector of labels, defining the classes to which the examples in
%in_features belong.
%
%in_k [default: 5]:
%k for k-nearest neighbors to be identified for each minority class index
%
%in_SMOTE [default: rounded ratio of majority class samples to minority class samples]
%positive integer, representing percentage of oversampling(assumption: integers = multiples of 100%, or 1). 
%

% OUTPUTS:
%----------
%out_featuresSyn, out_labelsSyn:
%features and labels of ONLY the synthetically created examples.
%note that each entry of out_labelsSyn is the label of the minority class
%since only examples of the minority class are created.
%concatenating [in_features out_featuresSyn] and [in_labels out_labelsSyn]
%gives a new example set with the desired class balance.

if nargin < 3 || isempty(in_k)
    in_k = 5;
end

mkNeg = 0; %handle -1s as 0s
if ~all(in_labels==0 | in_labels==1)
    if any(in_labels==-1)
        in_labels(in_labels==-1) = 0;
        mkNeg = 1;
    else
        error('SMOTE: in_labels may contain only the values 0 and 1.');
    end
end

numZeros = sum(in_labels==0);
numOnes  = sum(in_labels==1);

if numOnes == numZeros
    %nothing needs to be done because if classes are already balanced, then
    %for any in_beta, this is already the desired result.
    out_featuresSyn = [];
    out_labelsSyn   = [];
    return;
else
    if numZeros > numOnes
        majLabel = logical(0);
        minLabel = logical(1);
    else
        majLabel = logical(1);
        minLabel = logical(0);
    end
end

%rename:
S = in_features;
clear in_features;

%feature sets by class:
Smin = S(in_labels==minLabel,:);
Smaj = S(in_labels==majLabel,:);

%handle boundary cases:
if size(Smin,1)==0
    warning('SMOTE: there were no examples of the minority class in the data. hence balancing is not possible. Returning empty matrices.');
    out_featuresSyn = [];
    out_labelsSyn   = [];
    return;
end

%Ratio of samples in majority class
G = round(size(Smaj,1)/size(Smin,1))-1;

if size(Smin,1)==1
    warning('SMOTE: there was only one example of the minority class in the data. Hence returning G copies of that single example for balancing.');
    out_featuresSyn = repmat(Smin, [G 1]);
    out_labelsSyn   = logical(minLabel * ones([G 1]));
    return;
end

if nargin < 4 || isempty(in_SMOTE)
    in_SMOTE = G;
end

%determine nearest neighbors:
idcs = nearestneighbour(Smin',Smin','NumberOfNeighbours', in_k+1);
%note: why k+1? because we search kNNs of Smin in Smin itself and hence all
%points in Smin have a trivial nearest neighbor in Smin with distance 0.
%but that neighbor is not interesting because it's the point from Smin
%itself. hence remove it:
idcs = idcs';
idcs = idcs(:,2:end);

%initialize output and writing target as an empty matrix
Ssyn = zeros([0 size(Smin,2)]);

%for every minority example xi...
for cei=1:size(Smin,1)  %cei: current example index
    
    %current minority example:
    xi = Smin(cei,:);
    
    %number of synthetic examples to be created from xi:
    gi = in_SMOTE;
    
    %allocate space for gi examples to be created from xi:
    xiSyn = zeros(gi, size(Smin,2));
    
    %...iterate over synthetic examples to be created from xi and
    %random partner from set of nearest neighbors:
    for csi=1:gi    %csi: current synthetic example index
        
        %get random partner example from nearest neighbors of xi:
        %neighbor index:
        nIdx = idcs(cei, randi(size(idcs,2)));
        %neighbor:
        xiHat = Smin(nIdx,:);
        
        %create synthetic example as according to eq. (1) in reference [2]:
        delta = rand(1);
        xSyn = xi + delta * (xiHat - xi);
        
        %write it to xiSyn:
        xiSyn(csi,:) = xSyn;
        
    end
    
    %append examples synthesized from xi to overall synthetic example set:
    Ssyn = [Ssyn; xiSyn];
    
end

out_featuresSyn = Ssyn;
out_labelsSyn = (minLabel * ones([size(out_featuresSyn,1) 1]));

if mkNeg
    out_labelsSyn(out_labelsSyn==0) = -1;
else
    out_labelsSyn = logical(out_labelsSyn);
end

end

function [idx, tri] = nearestneighbour(varargin)
%NEARESTNEIGHBOUR    find nearest neighbours
%   IDX = NEARESTNEIGHBOUR(X) finds the nearest neighbour by Euclidean
%   distance to each point (column) in X from X. X is a matrix with points
%   as columns. IDX is a vector of indices into X, such that X(:, IDX) are
%   the nearest neighbours to X. e.g. the nearest neighbour to X(:, 2) is
%   X(:, IDX(2))
%
%   IDX = NEARESTNEIGHBOUR(P, X) finds the nearest neighbour by Euclidean
%   distance to each point in P from X. P and X are both matrices with the
%   same number of rows, and points are the columns of the matrices. Output
%   is a vector of indices into X such that X(:, IDX) are the nearest
%   neighbours to P
%
%   IDX = NEARESTNEIGHBOUR(I, X) where I is a logical vector or vector of
%   indices, and X has at least two rows, finds the nearest neighbour in X
%   to each of the points X(:, I).
%   I must be a row vector to distinguish it from a single point.
%   If X has only one row, the first input is treated as a set of 1D points
%   rather than a vector of indices
%
%   IDX = NEARESTNEIGHBOUR(..., Property, Value)
%   Calls NEARESTNEIGHBOUR with the indicated parameters set. Property
%   names can be supplied as just the first letters of the property name if
%   this is unambiguous, e.g. NEARESTNEIGHBOUR(..., 'num', 5) is equivalent
%   to NEARESTNEIGHBOUR(..., 'NumberOfNeighbours', 5). Properties are case
%   insensitive, and are as follows:
%      Property:                         Value:
%      ---------                         ------
%         NumberOfNeighbours             natural number, default 1
%            NEARESTNEIGHBOUR(..., 'NumberOfNeighbours', K) finds the closest
%            K points in ascending order to each point, rather than the
%            closest point. If Radius is specified and there are not
%            sufficient numbers, fewer than K neighbours may be returned
%
%         Radius                         positive, default +inf
%            NEARESTNEIGHBOUR(..., 'Radius', R) finds neighbours within
%            radius R. If NumberOfNeighbours is not set, it will find all
%            neighbours within R, otherwise it will find at most
%            NumberOfNeighbours. The IDX matrix is padded with zeros if not
%            all points have the same number of neighbours returned. Note
%            that specifying a radius means that the Delaunay method will
%            not be used.
%
%         DelaunayMode                   {'on', 'off', |'auto'|}
%            DelaunayMode being set to 'on' means NEARESTNEIGHBOUR uses the
%            a Delaunay triangulation with dsearchn to find the points, if
%            possible. Setting it to 'auto' means NEARESTNEIGHBOUR decides
%            whether to use the triangulation, based on efficiency. Note
%            that the Delaunay triangulation will not be used if a radius
%            is specified.
%
%         Triangulation                  Valid triangulation produced by
%                                        delaunay or delaunayn
%            If a triangulation is supplied, NEARESTNEIGHBOUR will attempt
%            to use it (in conjunction with dsearchn) to find the
%            neighbours.
%
%   [IDX, TRI] = NEARESTNEIGHBOUR( ... )
%   If the Delaunay Triangulation is used, TRI is the triangulation of X'.
%   Otherwise, TRI is an empty matrix
%
%   Example:
%
%     % Find the nearest neighbour in X to each column of X
%     x = rand(2, 10);
%     idx = nearestneighbour(x);
%
%     % Find the nearest neighbours to each point in p
%     p = rand(2, 5);
%     x = rand(2, 20);
%     idx = nearestneighbour(p, x)
%
%     % Find the five nearest neighbours to points x(:, [1 6 20]) in x
%     x = rand(4, 1000)
%     idx = nearestneighbour([1 6 20], x, 'NumberOfNeighbours', 5)
%
%     % Find all neighbours within radius of 0.1 of the points in p
%     p = rand(2, 10);
%     x = rand(2, 100);
%     idx = nearestneighbour(p, x, 'r', 0.1)
%
%     % Find at most 10 nearest neighbours to point p from x within a
%     % radius of 0.2
%     p = rand(1, 2);
%     x = rand(2, 30);
%     idx = nearestneighbour(p, x, 'n', 10, 'r', 0.2)
%
%
%   See also DELAUNAYN, DSEARCHN, TSEARCH

%TODO    Allow other metrics than Euclidean distance
%TODO    Implement the Delaunay mode for multiple neighbours

% Copyright 2006 Richard Brown. This code may be freely used and
% distributed, so long as it maintains this copyright line
narginchk(1, Inf);

% Default parameters
userParams.NumberOfNeighbours = []    ; % Finds one
userParams.DelaunayMode       = 'auto'; % {'on', 'off', |'auto'|}
userParams.Triangulation      = []    ;
userParams.Radius             = inf   ;

% Parse inputs
[P, X, fIndexed, userParams] = parseinputs(userParams, varargin{:});

% Special case uses Delaunay triangulation for speed.

% Determine whether to use Delaunay - set fDelaunay true or false
nX  = size(X, 2);
nP  = size(P, 2);
dim = size(X, 1);

switch lower(userParams.DelaunayMode)
    case 'on'
        %TODO Delaunay can't currently be used for finding more than one
        %neighbour
        fDelaunay = userParams.NumberOfNeighbours == 1 && ...
            size(X, 2) > size(X, 1)                    && ...
            ~fIndexed                                  && ...
            userParams.Radius == inf;
    case 'off'
        fDelaunay = false;
    case 'auto'
        fDelaunay = userParams.NumberOfNeighbours == 1 && ...
            ~fIndexed                                  && ...
            size(X, 2) > size(X, 1)                    && ...
            userParams.Radius == inf                   && ...
            ( ~isempty(userParams.Triangulation) || delaunaytest(nX, nP, dim) );
end

% Try doing Delaunay, if fDelaunay.
fDone = false;
if fDelaunay
    tri = userParams.Triangulation;
    if isempty(tri)
        try
            tri   = delaunayn(X');
        catch
            msgId = 'NearestNeighbour:DelaunayFail';
            msg = ['Unable to compute delaunay triangulation, not using it. ',...
                'Set the DelaunayMode parameter to ''off'''];
            warning(msgId, msg);
        end
    end
    if ~isempty(tri)
        try
            idx = dsearchn(X', tri, P')';
            fDone = true;
        catch
            warning('NearestNeighbour:DSearchFail', ...
                'dsearchn failed on triangulation, not using Delaunay');
        end
    end
else % if fDelaunay
    tri = [];
end

% If it didn't use Delaunay triangulation, find the neighbours directly by
% finding minimum distances
if ~fDone
    idx = zeros(userParams.NumberOfNeighbours, size(P, 2));

    % Loop through the set of points P, finding the neighbours
    Y = zeros(size(X));
    for iPoint = 1:size(P, 2)
        x = P(:, iPoint);

        % This is the faster than using repmat based techniques such as
        % Y = X - repmat(x, 1, size(X, 2))
        for i = 1:size(Y, 1)
            Y(i, :) = X(i, :) - x(i);
        end

        % Find the closest points, and remove matches beneath a radius
        dSq = sum(abs(Y).^2, 1);
        iRad = find(dSq < userParams.Radius^2);
        if ~fIndexed
            iSorted = iRad(minn(dSq(iRad), userParams.NumberOfNeighbours));
        else
            iSorted = iRad(minn(dSq(iRad), userParams.NumberOfNeighbours + 1));
            iSorted = iSorted(2:end);
        end

        % Remove any bad ones
        idx(1:length(iSorted), iPoint) = iSorted';
    end
    %while ~isempty(idx) && isequal(idx(end, :), zeros(1, size(idx, 2)))
    %    idx(end, :) = [];
    %end
    idx( all(idx == 0, 2), :) = [];
end % if ~fDone
if isvector(idx)
    idx = idx(:)';
end
end % nearestneighbour

%DELAUNAYTEST   Work out whether the combination of dimensions makes
%fastest to use a Delaunay triangulation in conjunction with dsearchn.
%These parameters have been determined empirically on a Pentium M 1.6G /
%WinXP / 512MB / Matlab R14SP3 platform. Their precision is not
%particularly important
function tf = delaunaytest(nx, np, dim)
switch dim
    case 2
        tf = np > min(1.5 * nx, 400);
    case 3
        tf = np > min(4 * nx  , 1200);
    case 4
        tf = np > min(40 * nx , 5000);

        % if the dimension is higher than 4, it is almost invariably better not
        % to try to use the Delaunay triangulation
    otherwise
        tf = false;
end % switch
end % delaunaytest

%MINN   find the n most negative elements in x, and return their indices
%  in ascending order
function I = minn(x, n)

% Make sure n is no larger than length(x)
n = min(n, length(x));

% Sort the first n
[xsn, I] = sort(x(1:n));

% Go through the rest of the entries, and insert them into the sorted block
% if they are negative enough
for i = (n+1):length(x)
    j = n;
    while j > 0 && x(i) < xsn(j)
        j = j - 1;
    end

    if j < n
        % x(i) should go into the (j+1) position
        xsn = [xsn(1:j), x(i), xsn((j+1):(n-1))];
        I   = [I(1:j), i, I((j+1):(n-1))];
    end
end

end %minn

%PARSEINPUTS    Support function for nearestneighbour
function [P, X, fIndexed, userParams] = parseinputs(userParams, varargin)
if length(varargin) == 1 || ~isnumeric(varargin{2})
    P           = varargin{1};
    X           = varargin{1};
    fIndexed    = true;
    varargin(1) = [];
else
    P             = varargin{1};
    X             = varargin{2};
    varargin(1:2) = [];

    % Check the dimensions of X and P
    if size(X, 1) ~= 1
        % Check to see whether P is in fact a vector of indices
        if size(P, 1) == 1
            try
                P = X(:, P);
            catch
                error('NearestNeighbour:InvalidIndexVector', ...
                    'Unable to index matrix using index vector');
            end
            fIndexed = true;
        else
            fIndexed = false;
        end % if size(P, 1) == 1
    else % if size(X, 1) ~= 1
        fIndexed = false;
    end

    if ~fIndexed && size(P, 1) ~= size(X, 1)
        error('NearestNeighbour:DimensionMismatch', ...
            'No. of rows of input arrays doesn''t match');
    end
end
% Parse the Property/Value pairs
if rem(length(varargin), 2) ~= 0
    error('NearestNeighbour:propertyValueNotPair', ...
        'Additional arguments must take the form of Property/Value pairs');
end

propertyNames = {'numberofneighbours', 'delaunaymode', 'triangulation', ...
    'radius'};
while length(varargin) ~= 0
    property = varargin{1};
    value    = varargin{2};

    % If the property has been supplied in a shortened form, lengthen it
    iProperty = find(strncmpi(property, propertyNames, length(property)));
    if isempty(iProperty)
        error('NearestNeighbour:InvalidProperty', 'Invalid Property');
    elseif length(iProperty) > 1
        error('NearestNeighbour:AmbiguousProperty', ...
            'Supplied shortened property name is ambiguous');
    end
    property = propertyNames{iProperty};

    switch property
        case 'numberofneighbours'
            if rem(value, 1) ~= 0 || ...
                    value > length(X) - double(fIndexed) || ...
                    value < 1
                error('NearestNeighbour:InvalidNumberOfNeighbours', ...
                    'Number of Neighbours must be an integer, and smaller than the no. of points in X');
            end
            userParams.NumberOfNeighbours = value;

        case 'delaunaymode'
            fOn = strcmpi(value, 'on');
            if strcmpi(value, 'off')
                userParams.DelaunayMode = 'off';
            elseif fOn || strcmpi(value, 'auto')
                if userParams.NumberOfNeighbours ~= 1
                    if fOn
                        warning('NearestNeighbour:TooMuchForDelaunay', ...
                            'Delaunay Triangulation method works only for one neighbour');
                    end
                    userParams.DelaunayMode = 'off';
                elseif size(X, 2) < size(X, 1) + 1
                    if fOn
                        warning('NearestNeighbour:TooFewDelaunayPoints', ...
                            'Insufficient points to compute Delaunay triangulation');
                    end
                    userParams.DelaunayMode = 'off';

                elseif size(X, 1) == 1
                    if fOn
                        warning('NearestNeighbour:DelaunayDimensionOne', ...
                            'Cannot compute Delaunay triangulation for 1D input');
                    end
                    userParams.DelaunayMode = 'off';
                else
                    userParams.DelaunayMode = value;
                end
            else
                warning('NearestNeighbour:InvalidOption', ...
                    'Invalid Option');
            end % if strcmpi(value, 'off')

        case 'radius'
            if isscalar(value) && isnumeric(value) && isreal(value) && value > 0
                userParams.Radius = value;
                if isempty(userParams.NumberOfNeighbours)
                    userParams.NumberOfNeighbours = size(X, 2) - double(fIndexed);
                end
            else
                error('NearestNeighbour:InvalidRadius', ...
                    'Radius must be a positive real number');
            end
    

        case 'triangulation'
            if isnumeric(value) && size(value, 2) == size(X, 1) + 1 && ...
                    all(ismember(1:size(X, 2), value))
                userParams.Triangulation = value;
            else
                error('NearestNeighbour:InvalidTriangulation', ...
                    'Triangulation not a valid Delaunay Triangulation');
            end
    end % switch property

    varargin(1:2) = [];
end % while
if isempty(userParams.NumberOfNeighbours)
    userParams.NumberOfNeighbours = 1;
end
end %parseinputs
