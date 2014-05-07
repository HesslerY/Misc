

function fr = fres(theta, phi, n, x, y, z, I)

fr = 0;
B = 2*pi;

for j = 1:n
    
    dr = -(x(j).*cos(phi).*sin(theta) + y(j).*sin(phi).*sin(theta) + z(j).*cos(theta));
    
    fr = (abs(I(j)).*exp(i .* angle(I(j))).*exp(-i .* B .* dr)) + fr;

end