% ============================================================================
%  add_lw_to_external.m  --  fill the LW column of rstar_external.csv from the
%  official Laubach-Williams (New York Fed) current-estimates workbook.
%
%  WHY:  Fig 6 overlays two standard r* measures on the DSGE forwards. Lubik-
%  Matthes (Richmond Fed) ships pre-filled in rstar_external.csv. LW is a binary
%  Excel file, so download it once and run this to add its column.
%
%  STEP 1 -- download the official file (updated ~quarterly):
%     https://www.newyorkfed.org/medialibrary/media/research/economists/williams/data/Laubach_Williams_current_estimates.xlsx
%     Save it next to this script.
%
%  STEP 2 -- from this folder:
%     >> add_lw_to_external
%     (or, if the auto-detect misses:  add_lw_to_external('Laubach_Williams_current_estimates.xlsx','United States') )
%
%  It writes the US natural rate (annual %) into rstar_external.csv, aligned by
%  quarter. Then re-run make_paper_figures to see LW on Fig 6.
%
%  The LW layout shifts between vintages, so this prints what it detected. If
%  the numbers look wrong, open the xlsx, note the sheet name + the column whose
%  header is 'rstar' (or 'r*') for the US, and pass them:
%     >> add_lw_to_external('...xlsx', SHEET, DATECOL, RSTARCOL)
% ============================================================================
function add_lw_to_external(XLSX, SHEET, DATECOL, RSTARCOL)
if nargin<1 || isempty(XLSX); XLSX = 'Laubach_Williams_current_estimates.xlsx'; end
assert(exist(XLSX,'file')==2, 'LW file not found: %s (download it first -- see header).', XLSX);
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
fprintf('LW: reading sheet "%s" from %s\n', char(SHEET), XLSX);
C = read_sheet_robust(XLSX, SHEET);

% ---- locate the header row + Date / rstar columns (auto), unless overridden ----
% The previous version took the FIRST cell matching /date/ and the first matching
% /rstar|natural rate/ and accepted them even when they were the SAME column --
% which is what happened on the LW workbook, whose row 1 is a title containing
% both words. It then read the date text as the value column and parsed nothing.
% This version requires dcol ~= rcol and VALIDATES the choice against the data
% below the header: the value column must be mostly finite numbers, and the date
% column must mostly parse to real quarters.
if nargin<4 || isempty(DATECOL) || isempty(RSTARCOL)
  [hdrRow,dcol,rcol] = find_cols(C);
  if hdrRow==0
    preview(C);
    error(['Could not auto-find usable Date/rstar columns. The sheet layout is ' ...
           'printed above: note the row holding the headers and the 1-based column ' ...
           'numbers, then call\n' ...
           '   add_lw_to_external(''%s'', ''%s'', DATECOL, RSTARCOL)'], XLSX, char(SHEET));
  end
  fprintf('LW: header row %d,  Date=col %d,  rstar=col %d\n', hdrRow, dcol, rcol);
else
  hdrRow = 1; dcol = DATECOL; rcol = RSTARCOL;
  fprintf('LW: using caller-supplied Date=col %d, rstar=col %d\n', dcol, rcol);
end

% ---- pull (date, rstar) pairs -> map to YYYYQq ------------------------------
vals = containers.Map('KeyType','char','ValueType','double');
badDates = {}; nNoNum = 0;
for i = hdrRow+1:size(C,1)
  dv = C{i,dcol}; rv = C{i,rcol};
  if ~(isnumeric(rv) && isscalar(rv) && isfinite(rv)); nNoNum = nNoNum+1; continue; end
  [Y,Q] = parse_lw_date(dv);
  if isnan(Y)
    % a row with a good NUMBER but an unparseable DATE is the dangerous case:
    % it drops silently and leaves a hole in the middle of the series.
    if ischar(dv)||isstring(dv); badDates{end+1} = char(dv); end %#ok<AGROW>
    continue
  end
  vals(sprintf('%dQ%d',Y,Q)) = rv;
end
fprintf('LW: parsed %d quarterly values\n', vals.Count);
if ~isempty(badDates)
  u = unique(badDates);
  fprintf('LW: *** %d rows had a valid number but an UNPARSEABLE DATE ***\n', numel(badDates));
  fprintf('LW: those become holes in the series. Distinct examples:\n');
  for k = 1:min(6,numel(u)); fprintf('        "%s"\n', u{k}); end
  % name the month tokens that failed, so the fix is obvious
  toks = regexp(badDates, '[A-Za-z\x80-\xFF]{3,}', 'match','once');
  toks = unique(lower(toks(~cellfun(@isempty,toks))));
  if ~isempty(toks)
    fprintf('LW: unrecognised month token(s): %s\n', strjoin(toks,', '));
    fprintf('LW: add them to the en_month map at the bottom of this file.\n');
  end
end
assert(vals.Count>0, 'No LW values parsed -- check the sheet/columns.');

% ---- rewrite rstar_external.csv, filling the LW column by date -------------
T = readcell(CSV);  hdr = T(1,:);
hcol = find(strcmpi(strtrim(string(hdr)),'LW_rstar'),1);
if isempty(hcol)   % accept a pre-rename CSV so an old copy is not given a 4th column
  hcol = find(strcmpi(strtrim(string(hdr)),'HLW_rstar'),1);
  if ~isempty(hcol); T{1,hcol} = 'LW_rstar'; end
end
if isempty(hcol); hcol = size(T,2)+1; T{1,hcol}='LW_rstar'; end
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
fprintf('LW: wrote LW column into %s. Re-run make_paper_figures for Fig 6.\n', CSV);
end

function [Y,Q] = parse_lw_date(dv)
% Returns NaN when the date cannot be read. It MUST NOT guess.
% The previous version fell back to "take the 4-digit year, call it Q1", which
% silently mapped every unparseable July row onto YYYYQ1 and OVERWROTE the real
% Q1 value: the series then lost Q3 and reported the Q3 number as Q1. Wrong data
% is far worse than missing data, and it hid from the parse-failure diagnostic
% because the fallback always returned a valid-looking year.
Y=NaN; Q=NaN;
if isnumeric(dv) && isscalar(dv) && isfinite(dv)
  if     dv>1900 && dv<2100;      Y=floor(dv+1e-6); Q=round((dv-Y)*4)+1;
  elseif dv>10000 && dv<100000;   d=datetime(dv,'ConvertFrom','excel');  Y=year(d); Q=quarter(d);
  elseif dv>700000;               d=datetime(dv,'ConvertFrom','datenum');Y=year(d); Q=quarter(d);
  end
elseif isdatetime(dv)
  Y=year(dv); Q=quarter(dv);
elseif ischar(dv)||isstring(dv)
  s=strtrim(char(dv));
  t=regexp(s,'(\d{4})\s*[:\-\.\s]?\s*[Qq]\s*([1-4])','tokens','once');      % 1961Q1
  if ~isempty(t); Y=str2double(t{1}); Q=str2double(t{2});
  else
    t=regexp(s,'^(\d{4})[-/](\d{1,2})[-/](\d{1,2})','tokens','once');       % ISO
    if ~isempty(t)
      Y=str2double(t{1}); Q=ceil(str2double(t{2})/3);
    else
      % dd-MON-yyyy in ANY language: pull the pieces out and map the month
      % ourselves rather than handing the string to datetime and hoping the
      % locale agrees.
      t=regexp(s,'(\d{1,2})\s*[-/\.\s]\s*([^\s\-/\.\d]+)\s*[-/\.\s]\s*(\d{2,4})','tokens','once');
      if ~isempty(t)
        M=month_num(t{2}); yy=str2double(t{3});
        if yy<100; if yy>=50; yy=yy+1900; else; yy=yy+2000; end; end
        if ~isnan(M); Y=yy; Q=ceil(M/3); end
      else
        d=[]; try; d=datetime(s); catch; end          % last resort, current locale
        if ~isempty(d); Y=year(d); Q=quarter(d); end
      end
    end
  end
end
if ~isnan(Q); Q=min(max(Q,1),4); end
end

function m = month_num(tok)
% Month number from an abbreviation in en/fr/it, accents and trailing dots
% removed first so encoding never matters: 'aout'/'aoUt' both reduce to 'aot',
% 'fevr'/'fEvr' to 'fvr', and so on. Longer keys are tested first because
% 'juin' and 'juil' share three letters, and 'juill' is a real French variant.
tok = lower(regexprep(char(tok),'[^A-Za-z]',''));
tab = { 'juill',7; 'juil',7; 'juin',6; 'janv',1; 'fevr',2; 'febr',2; 'mars',3; ...
        'aout',8; 'sept',9; ...
        'jan',1;  'feb',2;  'fvr',2;  'mar',3;  'avr',4;  'apr',4;  'mai',5; ...
        'may',5;  'mag',5;  'jun',6;  'giu',6;  'jul',7;  'lug',7;  'aot',8; ...
        'aug',8;  'ago',8;  'sep',9;  'set',9;  'oct',10; 'ott',10; 'nov',11; ...
        'dec',12; 'dic',12; 'dc',12 };
m = NaN;
for k = 1:size(tab,1)
  pk = tab{k,1};
  if numel(tok)>=numel(pk) && strncmp(tok,pk,numel(pk)); m = tab{k,2}; return; end
end
end

function C = read_sheet_robust(XLSX, SHEET)
% readcell converts date-looking cells to datetime using the CURRENT locale and
% throws on e.g. '01-janv-1961' when the locale is en_US. Ask for raw Excel
% serial numbers instead: locale-proof, and parse_lw_date handles serials.
  C = [];
  try
    C = readcell(XLSX,'Sheet',SHEET,'DatetimeType','exceldatenum');
    fprintf('LW: read with DatetimeType=exceldatenum (locale-proof)\n'); return
  catch
  end
  try
    C = readcell(XLSX,'Sheet',SHEET,'DatetimeType','text');
    fprintf('LW: read with DatetimeType=text\n'); return
  catch
  end
  C = readcell(XLSX,'Sheet',SHEET);
  fprintf('LW: read with default settings\n');
end

