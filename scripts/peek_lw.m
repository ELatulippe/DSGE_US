% ============================================================================
%  peek_lw.m  --  show what is actually in the LW workbook, column by column.
%
%  add_lw_to_external picked column 3 and got 197 values with NO date-parsing
%  failures, which means the 64 missing rows simply had no number in that
%  column. So column 3 is very likely not the series we want -- the LW workbook
%  carries several (one-sided, two-sided, trend growth g, the z component), and
%  they do not all span the same rows.
%
%  This prints, for every column: its header, how many finite numbers it holds,
%  the first and last quarter it covers, and its value range. Pick the column
%  whose count is largest, whose range ends in the current year, and whose
%  values sit where r* should (LW r* is roughly 0.5 to 3.5 percent), then:
%
%     >> add_lw_to_external('Laubach_Williams_current_estimates.xlsx','data',1,COL)
%
%     >> peek_lw
%     >> peek_lw('Laubach_Williams_current_estimates.xlsx','data',6)
% ============================================================================
function peek_lw(XLSX, SHEET, HDRROW)

if nargin<1 || isempty(XLSX);   XLSX  = 'Laubach_Williams_current_estimates.xlsx'; end
if nargin<2 || isempty(SHEET);  SHEET = 'data'; end
if nargin<3 || isempty(HDRROW); HDRROW = 6; end
assert(exist(XLSX,'file')==2,'%s not found in this folder.',XLSX);

C = [];
for dt = {'exceldatenum','text'}
  try; C = readcell(XLSX,'Sheet',SHEET,'DatetimeType',dt{1}); break; catch; end
end
if isempty(C); C = readcell(XLSX,'Sheet',SHEET); end
[nR,nC] = size(C);
fprintf('\n===== %s / sheet "%s" : %d rows x %d cols, header row %d =====\n\n', ...
        XLSX, SHEET, nR, nC, HDRROW);

% ---- the header rows, so you can see the layout ----------------------------
fprintf('--- rows 1..%d ---\n', min(HDRROW+2,nR));
for i = 1:min(HDRROW+2,nR)
  fprintf('r%-3d|', i);
  for j = 1:min(12,nC); fprintf(' %-14s', cellstr_(C{i,j})); end
  fprintf('\n');
end

% ---- per-column census -----------------------------------------------------
fprintf('\n--- per-column census (rows %d..%d) ---\n', HDRROW+1, nR);
fprintf('%4s  %-26s %7s  %-9s %-9s %10s %10s\n', ...
        'col','header','numbers','first','last','min','max');
fprintf('%s\n', repmat('-',1,88));
for j = 1:nC
  hdr = cellstr_(C{HDRROW,j});
  v = nan(nR,1); q = cell(nR,1);
  for i = HDRROW+1:nR
    x = C{i,j};
    if isnumeric(x) && isscalar(x) && isfinite(x); v(i) = x; end
    [Y,Q] = pdate(C{i,1});
    if ~isnan(Y); q{i} = sprintf('%dQ%d',Y,Q); end
  end
  ok = ~isnan(v);
  if ~any(ok); continue; end
  qq = q(ok); qq = qq(~cellfun(@isempty,qq));
  f = ''; l = '';
  if ~isempty(qq); f = qq{1}; l = qq{end}; end
  fprintf('%4d  %-26s %7d  %-9s %-9s %10.4f %10.4f\n', ...
          j, trunc(hdr,26), sum(ok), f, l, min(v(ok)), max(v(ok)));
end
fprintf('%s\n', repmat('-',1,88));
fprintf('\nPick the column with the MOST numbers, a last quarter in the current year,\n');
fprintf('and values in a plausible r* range (LW r* runs about 0.5 to 3.5 percent).\n');
fprintf('Then:  add_lw_to_external(''%s'', ''%s'', 1, COL)\n\n', XLSX, char(SHEET));
end

% ---- helpers ---------------------------------------------------------------
function s = cellstr_(x)
  if ischar(x)||isstring(x); s = char(x);
  elseif isnumeric(x)&&isscalar(x)&&isfinite(x); s = sprintf('%.6g',x);
  else; s = '.'; end
  s = trunc(s,14);
end
function s = trunc(s,n)
  s = strtrim(s); if numel(s)>n; s = [s(1:n-1) '~']; end
end
function [Y,Q] = pdate(dv)
Y=NaN; Q=NaN;
if isnumeric(dv) && isscalar(dv) && isfinite(dv)
  if     dv>1900 && dv<2100;      Y=floor(dv+1e-6); Q=round((dv-Y)*4)+1;
  elseif dv>10000 && dv<100000;   d=datetime(dv,'ConvertFrom','excel'); Y=year(d); Q=quarter(d);
  elseif dv>700000;               d=datetime(dv,'ConvertFrom','datenum'); Y=year(d); Q=quarter(d);
  end
elseif isdatetime(dv); Y=year(dv); Q=quarter(dv);
elseif ischar(dv)||isstring(dv)
  s=char(dv);
  t=regexp(s,'(\d{4})\s*[:\-\.\s]?\s*[Qq]\s*([1-4])','tokens','once');
  if ~isempty(t); Y=str2double(t{1}); Q=str2double(t{2}); return; end
  d=[]; try; d=datetime(s); catch; end
  if ~isempty(d); Y=year(d); Q=quarter(d);
  else; n=regexp(s,'\d{4}','match','once'); if ~isempty(n); Y=str2double(n); Q=1; end; end
end
if ~isnan(Q); Q=min(max(Q,1),4); end
end
