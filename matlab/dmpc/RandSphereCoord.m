function pos = RandSphereCoord(r,c,zlim)
    % rand(1)*2-1 gives -1 to 1
    theta = (rand(1)*2-1)*pi;
    phi = (rand(1)*2*zlim-zlim)*pi/2;
    cosphi = cos(phi);
    sintheta = sin(theta);
    pos(1) = r(1)*cosphi*cos(theta)  + c(1);
    pos(2) = r(2)*cosphi*sintheta  + c(2);
    pos(3) = r(3)*sin(phi)  + c(3);
end