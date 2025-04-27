function [axpos] = subplotpos_col(nrow, ncol, varargin)
% SUBPLOTPOS_COL returns axes positions for nrow x ncol subplots.
% Populates columns first.
narginchk(2, 10);
validateattributes(nrow, {'numeric'}, {'scalar', 'integer', 'positive'});
validateattributes(ncol, {'numeric'}, {'scalar', 'integer', 'positive'});
npnl = nrow * ncol;
p = inputParser;
p.addParameter('xgap', 0.002, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', '<=', 0.5}));
p.addParameter('ygap', 0.002, @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', '<=', 0.5}));
p.addParameter('lr', [0.08 0.92], @(x) validateattributes(x, {'numeric'}, {'increasing', 'numel', 2, '>=', 0, '<=', 1}));
p.addParameter('bu', [0.079 0.94], @(x) validateattributes(x, {'numeric'}, {'increasing', 'numel', 2, '>=', 0, '<=', 1}));
parse(p, varargin{:});
xgap = p.Results.xgap;
ygap = p.Results.ygap;
lr = p.Results.lr;
bu = p.Results.bu;
fig_w = lr(2) - lr(1);
fig_h = bu(2) - bu(1);
pw = (1 - xgap) * fig_w / ncol;
ph = (1 - ygap) * fig_h / nrow;
mx = linspace(lr(1), lr(2), 2 * ncol + 1);
mx = mx(2:2:end);
my = linspace(bu(1), bu(2), 2 * nrow + 1);
my = my(2:2:end);
xo = repmat(mx - pw / 2, nrow, 1);
yo = flipud(repmat((my - ph / 2).', 1, ncol));
wz = pw * ones(npnl, 1);
hz = ph * ones(npnl, 1);
axpos = [xo(:), yo(:), wz, hz];
end
