b = [0, 0, 0, 1];
a = [1, 12, 40, 0];

H = tf(b, a);
tvals = linspace(0, 10, 100);

u = 1 - (sqrt(10)/4)*exp(-6*tvals).*sin(2*tvals + atan(1/3));

figure; 
step(H)