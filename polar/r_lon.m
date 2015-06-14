% copy the final computed kernel to this folder. Do the ncessary changes in this code according to your requirement. Run this and get the plot
% it plots radius v/s longitude


f = load('JCD_modelfull');
r = f(:,1);
r=flipud(r);
r= r(1:8:1944);

load 'Dkernel.mat';

for ir = 75:272
    for ith = 125:175

        rad = r(ir);
        th = (ith-1)*(pi/300);

        x(ir-74,ith) = rad * cos(th);
        y(ir-74,ith) = rad * sin(th);
        c(ir-74,ith) = Dkernel(ir,150,ith);
    end
end
%x(:,302)=[];
%y(:,302)=[];
size(x)
size(y)
size(c)
hold on
pcolor(x,y,c)
shading interp
%hold on
%polarmap(x,y,c,'k-')
colorbar
axis equal tight

