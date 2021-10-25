function [sources_int, sources_nonint, P_ar] = generate_sources_ar(fs, len, bandpass)
% Stefan Haufe, 2014, 2015
% stefan.haufe@tu-berlin.de
%
% If you use this code for a publication, please ask Stefan Haufe for the
% correct reference to cite.

% License
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see http://www.gnu.org/licenses/.

N = fs*len;
P_ar = 5;

rep = 1;
while rep
  sources_nonint = zeros(N+100, 2);
  sources_nonint(:, 1) = gen_ar_uni(N+100, P_ar);
  sources_nonint(:, 2) = gen_ar_uni(N+100, P_ar);
  if ~isempty(bandpass)   
    [P1, f1] = pwelch(sources_nonint(:, 1), hanning(fs), [], fs, fs);
    [P2, f2] = pwelch(sources_nonint(:, 2), hanning(fs), [], fs, fs);
    P = P1 + P2;
    in_f = find(ismember(f1, bandpass(1):bandpass(2)));
    if (sum(P(in_f))/sum(P) > 1.2*length(in_f)/length(f1))  
      [b, a] = butter(3, bandpass/fs*2);
      sources_nonint = filtfilt(b, a, sources_nonint);
      rep = 0;
    end
  else
    rep = 0;      
  end
  sources_nonint = sources_nonint((100+1):end, :)';
end

rep = 1;
while rep
  sources_int = gen_ar_biv(N+100, P_ar)';

  if ~isempty(bandpass)   
    [P1, f1] = pwelch(sources_int(:, 1), hanning(fs), [], fs, fs);
    [P2, f2] = pwelch(sources_int(:, 2), hanning(fs), [], fs, fs);
    P = P1 + P2;
    in_f = find(ismember(f1, bandpass(1):bandpass(2)));
    if (sum(P(in_f))/sum(P) > 1.2*length(in_f)/length(f1))  
      [b, a] = butter(3, bandpass/fs*2);
      sources_int = filtfilt(b, a, sources_int);
      rep = 0;
    end
  else
    rep = 0;      
  end
  sources_int = sources_int((100+1):end, :)';
end

end

