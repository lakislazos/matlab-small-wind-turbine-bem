function polar = loadPolar(csvFile)
% LOADPOLAR  Read an aerofoil polar CSV into a struct with interpolants.
%
%   polar = loadPolar(csvFile)
%
% The CSV is expected to have columns named (case-insensitive):
%   alpha (deg), CL, CD     and optionally CM
%
% This handles the standard AirfoilTools export format, which has a
% multi-line header followed by the data. We try to detect the header
% row automatically by finding the first row containing the word 'alpha'.
%
% RETURNS a struct with:
%   polar.alpha      : column vector of angles of attack (deg)
%   polar.CL         : column vector of lift coefficients
%   polar.CD         : column vector of drag coefficients
%   polar.CL_fn      : griddedInterpolant for C_L(alpha)
%   polar.CD_fn      : griddedInterpolant for C_D(alpha)
%   polar.alphaMin   : min alpha in data (for clamping)
%   polar.alphaMax   : max alpha in data
%   polar.file       : original filename (for logging)

    % --- Find the header line ---
    fid = fopen(csvFile, 'r');
    if fid < 0
        error('loadPolar:fileNotFound', 'Cannot open %s', csvFile);
    end
    cleanup = onCleanup(@() fclose(fid));

    headerRow = -1;
    lineIdx = 0;
    while ~feof(fid)
        line = fgetl(fid);
        lineIdx = lineIdx + 1;
        if ischar(line) && contains(lower(line), 'alpha')
            headerRow = lineIdx;
            break;
        end
    end
    if headerRow < 0
        error('loadPolar:noHeader', ...
              'Could not find a header row containing "alpha" in %s', csvFile);
    end

    % --- Read the table starting at the header line ---
    opts = detectImportOptions(csvFile, 'NumHeaderLines', headerRow - 1);
    opts.VariableNamingRule = 'preserve';
    T = readtable(csvFile, opts);

    % Normalize column names: lower-case, strip whitespace
    cols = lower(strtrim(T.Properties.VariableNames));
    T.Properties.VariableNames = cols;

    % Find the columns we need
    alphaCol = findCol(cols, {'alpha', 'aoa', 'angle'});
    clCol    = findCol(cols, {'cl', 'c_l'});
    cdCol    = findCol(cols, {'cd', 'c_d'});

    alpha = T.(alphaCol);
    CL    = T.(clCol);
    CD    = T.(cdCol);

    % Drop any rows with NaN
    keep = ~(isnan(alpha) | isnan(CL) | isnan(CD));
    alpha = alpha(keep);
    CL    = CL(keep);
    CD    = CD(keep);

    % Sort by alpha (required for griddedInterpolant)
    [alpha, idx] = sort(alpha);
    CL = CL(idx);
    CD = CD(idx);

    % Remove duplicates (some airfoil databases include alpha both
    % sweeping up and sweeping down)
    [alpha, ia] = unique(alpha, 'stable');
    CL = CL(ia);
    CD = CD(ia);

    % Build interpolants. 'linear' interpolation in-range, 'nearest'
    % extrapolation outside the data — this is a deliberate choice:
    % BEM iterations can briefly probe alpha values outside the polar
    % range during convergence, and using 'nearest' keeps the solver
    % stable rather than letting C_L blow up linearly.
    polar.CL_fn = griddedInterpolant(alpha, CL, 'linear', 'nearest');
    polar.CD_fn = griddedInterpolant(alpha, CD, 'linear', 'nearest');

    polar.alpha    = alpha;
    polar.CL       = CL;
    polar.CD       = CD;
    polar.alphaMin = min(alpha);
    polar.alphaMax = max(alpha);
    polar.file     = csvFile;
end


function colName = findCol(cols, candidates)
% Find the first column whose name matches one of the candidates.
    for k = 1:numel(candidates)
        idx = find(strcmp(cols, candidates{k}), 1);
        if ~isempty(idx)
            colName = cols{idx};
            return;
        end
    end
    error('loadPolar:missingCol', ...
          'Could not find any of [%s] in columns [%s]', ...
          strjoin(candidates, ', '), strjoin(cols, ', '));
end
