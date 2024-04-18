function res = mph2fpm(xInp, j)
  c = 1.605/3.6/0.3048*60;
  if j > 0 
    res = xInp*c;
  else
    res = xInp/c;
  end
      