function yNew = rk4_Wadehn(t,y,f,h)
%RK4_WADEHN Summary of this function goes here
%   Detailed explanation goes here

k1 = h*f(t, y);
k2 = h*f(t+h/2, y+k1/2);
k3 = h*f(t+h/2, y+k2/2);
k4 = h*f(t+h, y+k3);

yNew = y + (k1+2*k2+2*k3+k4)/6;

end

