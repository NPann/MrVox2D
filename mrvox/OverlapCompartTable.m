% Copyright 2013 Nicolas Pannetier, Clément Debacker, Franck Mauconduit, Thomas Christen, Emmanuel Barbier
% This file is part of DCESIM.
% 
%   DCESIM is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%    DCESIM is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with DCESIM.  If not, see <http://www.gnu.org/licenses/>.

function OverlapMatrice = OverlapCompartTable(Centers,Radius)


C = Centers;
R = Radius;
n = size(C,1);

dist = ...
    (( C(:,1)*ones(1,n) - ones(n,1)*C(:,1)' ) .^2 + ...
     ( C(:,2)*ones(1,n) - ones(n,1)*C(:,2)' ) .^2 ).^(1/2);
 
dist(dist==0) = Inf;

OverlapMatrice = dist < ones(n,1)*R(1:n) + (ones(n,1)*R(1:n))';
