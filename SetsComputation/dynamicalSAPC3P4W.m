function dx = dynamicalSAPC3P4W(t,x0,u)

%Parametros del sistema

f = 50;
w = 2*pi*f;

Lp = 30.2e-3;
Rp = 1;

Cdc = 2200e-6;
Ro = 1000;

% Estado estacionario

ipd = x0(1);
ipq = x0(2);

vdc = x0(3);

upd = u(1);
upq = u(2);
vpccd = u(3);


dipd = (1/Lp)*(Lp*w*ipq-Rp*ipd-upd*(vdc/2)+vpccd);

dipq = (1/Lp)*(-Lp*w*ipd-Rp*ipq-upq*(vdc/2));

dvdc = (1/Cdc)*(upd*ipd+upq*ipq-(vdc/Ro));

dx = [dipd, dipq, dvdc]';

end

