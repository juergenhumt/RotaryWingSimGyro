 function [frcG] = frcGFunc( phi, thta)
%
% Copyright 2010/2011 Juergen Humt
% 
% This file is part of RotaryWingSim.
% 
%     RotaryWingSim, is free  software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by the 
%     Free Software Foundation, either version 3 of the License or any later 
%     version.
% 
%     RotaryWingSim is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License along 
%     with RotaryWingSim.  If not, see <http://www.gnu.org/licenses/>.
%
%
%  Version 1.0 of this program was based on report naca-tm-73254 
%
%         --- Version 1.2 ---
%
%
% This modul transforms the gravity force to body coordinates.
%
% Input:
% phi  : roll angle
% thta : pitch angle
%
% Output
% frcG : gravity force in body coordinates
% 
%

 global  mHeliG
 
  %
  fx=0; fy=0; fz=mHeliG;
  psi=0;
  
  frcG = RTrfFunc(phi, thta, psi, fx,fy,fz);
  frcG = frcG';