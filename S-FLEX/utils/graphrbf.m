function [D B] = graphrbf(surf, sigma, inds)
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


Nvc = size(surf.vc, 1);
Ntri = size(surf.tri, 1);
% G = sparse(surf.tri(:,[3 1 2]), surf.tri, 1, Nvc, Nvc);

GW = spalloc(Nvc, Nvc, 3*Ntri);
for ii = 1:Ntri
    tri = surf.tri(ii, :);
    GW(tri(1), tri(2)) = norm(surf.vc(tri(1), :)-surf.vc(tri(2), :));
    GW(tri(2), tri(3)) = norm(surf.vc(tri(2), :)-surf.vc(tri(3), :));
    GW(tri(3), tri(1)) = norm(surf.vc(tri(3), :)-surf.vc(tri(1), :)); 
end

if nargin > 2
  for iind = 1:length(inds)
    D(:, iind) = shortest_paths(GW, inds(iind));
  end
else
  D = johnson_all_sp(GW);
end

% D = graphallshortestpaths(GW);

if nargin > 1
    B = normpdf(D, 0, sigma);
else
    B = [];
end

end

