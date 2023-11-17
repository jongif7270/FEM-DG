function n = normal(v)
x=v(1);y=v(2);
R = [0 1 ; -1 0];
n = R*[x;y]/(x^2+y^2)^(1/2);
