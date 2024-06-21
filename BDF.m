function [w,gacorr] = BDF_A(f,T,a,b,d)

len = length(f);
fext = [zeros(1,length(f)),f,zeros(1,length(f))];

for i=1:1:len % Caculate the instantaneous auto-correlation function (IAF)
    for j = -(len-1)/2:1:(len-1)/2
        acorr(i,j+(len-1)/2+1) = fext(i+len+j)*conj(fext(i+len-j));
    end
end

gacorr=acorr';
w=gacorr;

for i=1:len % Calculate the lct of each line in the IAF 
    w(:,i) = direct_method1d(gacorr(:,i)', 2*T, a, b, d);
end
end