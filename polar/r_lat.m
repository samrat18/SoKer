% copy the final computed kernel to this folder. Do the ncessary changes in this code according to your requirement. Run this and get the plot

%it plot radius v/s latitude

f = load('JCD_modelfull');
r = f(:,1);
r=flipud(r);
r= r(1:8:1944);

load 'Dkernel.mat';

for ir = 75:243
    for iph = 125:225
       
        rad = r(ir); 
        ph = (iph-1)*(2*pi/600);

        x(ir-74,iph) = rad * cos(ph);
        y(ir-74,iph) = rad * sin(ph);
        c(ir-74,iph) = Dkernel(ir,iph,151);
    end
end
%x(:,302)=[];
%y(:,302)=[];

hold on
pcolor(x,y,c)
shading interp
%hold on
%polarmap(x,y,c,'k-')
colorbar
axis equal tight
