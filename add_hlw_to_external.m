% ============================================================================
%  add_hlw_to_external.m  --  fill the HLW column of rstar_external.csv from the
%  official Holston-Laubach-Williams (New York Fed) current-estimates workbook.
%
%  WHY:  Fig 6 overlays two standard r* measures on the DSGE forwards. Lubik-
%  Matthes (Richmond Fed) ships pre-filled in rstar_external.csv. HLW is a binary
%  Excel file, so download it once and run this to add its column.
%
%  STEP 1 -- download the official file (updated ~quarterly):
%     https://www.newyorkfed.org/medialibrary/media/research/economists/williams/data/Holston_Laubach_Williams_current_estimates.xlsx
%     Save it next to this script.
%
%  STEP 2 -- from this folder:
%     >> add_hlw_to_external
%     (or, if the auto-detect misses:  add_hlw_to_external('Holston_Laubach_Williams_current_estimates.xlsx','United States') )
%
%  It writes the US natural rate (annual %) into rstar_external.csv, aligned by
%  quarter. Then re-run make_paper_figures to see HLW on Fig 6.
%
%  The HLW layout shifts between vintages, so this prints what it detected. If
%  the numbers look wrong, open the xlsx, note the sheet name + the column whose
%  header is 'rstar' (or 'r*') for the US, and pass them:
%     >> add_hlw_to_external('...xlsx', SHEET, DATECOL, RSTARCOL)
% ============================================================================
function add_hlw_to_external(XLSX, SHEET, DATECOL, RSTARCOL)
if nargin<1 || isempty(XLSX); XLSX = 'Holston_Laubach_Williams_current_estimates.xlsx'; end
assert(exist(XLSX,'file')==2, 'HLW file not found: %s (download it first -- see header).', XLSX);
CSV = 'rstar_external.csv';
assert(exist(CSV,'file')==2, 'rstar_external.csv not found in this folder.');

% ---- pick the sheet (default: the one whose name mentions the United States) --
sheets = sheetnames(XLSX);
if nargin<2 || isempty(SHEET)
  hit = find(contains(lower(sheets),'united states') | contains(lower(sheets),'u.s') | strcmpi(sheets,'us'),1);
  if isempty(hit); hit = find(contains(lower(sheets),'data'),1); end
  if isempty(hit); hit = 1; end
  SHEET = sheets{hit};
end
fprintf('HLW: reading sheet "%s" from %s\n', char(SHEET), XLSX);
C = readcell(XLSX,'Sheet',SHEET);

% ---- locate the header row + Date / rstar columns (auto), unless overridden ----
if nargin<4 || isempty(DATECOL) || isempty(RSTARCOL)
  hdrRow = 0; dcol = 0; rcol = 0;
  for i = 1:min(12,size(C,1))
    row = C(i,:);
    isTxt = cellfun(@(x) ischar(x)||isstring(x), row);
    lc = repmat({''},1,numel(row)); lc(isTxt) = lower(strtrim(cellstr(string(row(isTxt)))));
    dc = find(strcmp(lc,'date') | contains(lc,'date'),1);
    rc = find(strcmp(lc,'rstar') | strcmp(lc,'r*') | contains(lc,'rstar') | contains(lc,'natural rate'),1);
    if ~isempty(dc) && ~isempty(rc); hdrRow=i; dcol=dc; rcol=rc; break; end
  end
  assert(hdrRow>0, 'Could not auto-find Date/rstar headers -- pass DATECOL, RSTARCOL explicitly.');
  fprintf('HLW: header row %d,  Date=col %d,  rstar=col %d\n', hdrRow, dcol, rcol);
else
  hdrRow = 1; dcol = DATECOL; rcol = RSTARCOL;
end

% ---- pull (date, rstar) pairs -> map to YYYYQq ------------------------------
vals = containers.Map('KeyType','char','ValueType','double');
for i = hdrRow+1:size(C,1)
  dv = C{i,dcol}; rv = C{i,rcol};
  if ~(isnumeric(rv) && isscalar(rv) && isfinite(rv)); continue; end
  [Y,Q] = parse_hlw_date(dv);
  if isnan(Y); continue; end
  vals(sprintf('%dQ%d',Y,Q)) = rv;
end
fprintf('HLW: parsed %d quarterly values\n', vals.Count);
assert(vals.Count>0, 'No HLW values parsed -- check the sheet/columns.');

% ---- rewrite rstar_external.csv, filling the HLW column by date -------------
T = readcell(CSV);  hdr = T(1,:);
hcol = find(strcmpi(strtrim(string(hdr)),'HLW_rstar'),1);
if isempty(hcol); hcol = size(T,2)+1; T{1,hcol}='HLW_rstar'; end
for i = 2:size(T,1)
  d = strtrim(char(string(T{i,1})));
  if isKey(vals,d); T{i,hcol} = vals(d); else; T{i,hcol} = []; end
end
fid = fopen(CSV,'w');
for i = 1:size(T,1)
  for j = 1:size(T,2)
    x = T{i,j};
    if j>1; fprintf(fid,','); end
    if ischar(x)||isstring(x); fprintf(fid,'%s',char(x));
    elseif isnumeric(x)&&isscalar(x)&&isfinite(x); fprintf(fid,'%.6f',x);
    end
  end
  fprintf(fid,'\n');
end
fclose(fid);
fprintf('HLW: wrote HLW column into %s. Re-run make_paper_figures for Fig 6.\n', CSV);
end

function [Y,Q] = parse_hlw_date(dv)
Y=NaN; Q=NaN;
if isnumeric(dv) && isscalar(dv)
  if dv>1900 && dv<2100                      % decimal year e.g. 2025.25
    Y=floor(dv+1e-6); Q=round((dv-Y)*4)+1;
  elseif dv>700000                            % Excel serial date
    dt=datetime(dv,'ConvertFrom','excel'); Y=year(dt); Q=quarter(dt);
  end
elseif isdatetime(dv)
  Y=year(dv); Q=quarter(dv);
elseif ischar(dv)||isstring(dv)
  s=char(dv); tok=regexp(s,'(\d{4}).*?([1-4Qq:\-\.]?)\s*Q?\s*([1-4])?','tokens','once');
  try; dt=datetime(s); Y=year(dt); Q=quarter(dt); catch
    n=regexp(s,'\d{4}','match','once'); if ~isempty(n); Y=str2double(n); Q=1; end
  end
end
if ~isnan(Q); Q=min(max(Q,1),4); end
end
