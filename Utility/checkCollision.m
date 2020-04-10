function isCollided = checkCollision(x1,y1,R1,x2,y2,R2)
[xout, yout] = circcirc(x1,y1,R1,x2,y2,R2);
if isnan(xout)
    isCollided = false;
else
    isCollided = true;
end