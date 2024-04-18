function printJcb(jcb,jDim,lgPrnt)
%
% Copyright 2010 Juergen Humt
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
%         --- Version 1.0 ---
%
%
% This modul prints the jacobian matrix of the
% aircraft trim modul
%
% Input:
% jcb   : jacobian matrix to print
% jDim  : dimension of jacobian matrix
% lgPrnt: switch to allow printing of various matrix header lines
% 
% Output
% none (screen only)
% 
%
%
  if 7==jDim
    if lgPrnt > 0
      if 1==lgPrnt  
        fprintf('             thtaN          A1           B1           thtaT         phiHeli       thtaHeli      oM\n');
      else
%   trmVec=[1 2 3 4 6 7 9]; 
%                         AcFrcMmtSummationOutF(dlc,      dle,      dla,      dlp,      dlt,     phi      thtaFus     oM       
% [fext mext oMDot fZMR]= AcFrcMmtSummationOutF(trmVr(1), trmVr(2), trmVr(3), trmVr(4), trmVr(5),trmVr(6), trmVr(7), trmVr(8), trmVr(9), vErth, omBd,  mHeliG, outP, outFile);
          
        fprintf('              u             v             w             p            q             r\n');
      end  
      charAr=['u.','v.','w.','p.','q.','r.'];      
    end
    if -1==lgPrnt
      fprintf('           dlc           dle           dla           dlp           phiFus        thtaFus       oM\n');
      charAr=['tN','A1','B1','tT','fi','tH','oM'];
    end  
  end  


  if 6==jDim
    if lgPrnt > 0
      if 1==lgPrnt  
        fprintf('             thtaN          A1           B1           thtaT         phiHeli       thtaHeli\n');
      else
        fprintf('              u             v             w             p            q             r\n');
      end  
      charAr=['u.','v.','w.','p.','q.','r.'];      
    end
    if -1==lgPrnt
      fprintf('               u.           v.           w.             p.            q.              r.\n');
      charAr=['tN','A1','B1','tT','fi','tH'];
    end  
  end  

if 5==jDim
    if 1==lgPrnt
      fprintf('             thtaN           A1            B1         thtaT       thtaHeli\n');
      charAr=['u.','w.','p.','q.','r.'];      
    end
    if -1==lgPrnt
      fprintf('               u.           w.            p.            q.           r.\n');
      charAr=['tN','A1','B1','tT','tH'];
    end  
  end  
  
  
  
  
 if 4==jDim
    if 1==lgPrnt
      fprintf('             thtaN           B1        thtaT       thtaHeli\n');
      charAr=['u.','w.','q.','r.'];      
    end
    if -1==lgPrnt
      fprintf('               u.           w.            q.           r.\n');
      charAr=['tN','B1','tT','tH'];
    end  
  end  

   
  if 3==jDim
    if lgPrnt > 0
      if 1==lgPrnt  
        fprintf('             thtaN           B1        thtaHeli\n');
      else
        fprintf('             u             w             q\n');
      end  
      charAr=['u.','w.','q.'];      
    end
    if -1==lgPrnt
      fprintf('               u.           w.             q.\n');
      charAr=['tN','B1','tH'];
    end  
  end  
  

  if 12==jDim
    if 1==lgPrnt
      fprintf('             thtaN             B1\n');
      charAr=['u.','w.'];      
    end
    if -1==lgPrnt
      fprintf('               u.           w.\n');
      charAr=['tN','B1'];
    end  
  end  
  
  
  if   2 ==jDim  % 2 !!!
    if 1==lgPrnt
      fprintf('             thtaN           thtaHeli\n');
      charAr=['u.','w.'];      
    end
    if -1==lgPrnt
      fprintf('               u.           w.\n');
      charAr=['tN','tH'];
    end  
  end
  
  if lgPrnt ~= 0 
  for i=1:jDim
     k=2*i;
     if jDim > 1
       fprintf('%d  %c%c ',i,charAr(k-1),charAr(k));
     end  
     for j=1:jDim-1;
         fprintf('  %12.3e',jcb(i,j));
     end
     j=j+1;
     fprintf('  %12.3e\n',jcb(i,j))
  end
  
  dt=det(jcb);
  fprintf('Det %12.7e\n',dt);
  end