% 

function [f, df] = sampleEijOpt(x)
% optimization code

% param
vert = x(4:7);
w_res = x(8);

% box case (4 lines), three free parameters
x0 = x(1); y0 = x(2); % camera center
yc = x(3); % corner

% assume vertical lines are fixed
a_aob = (vert(2) - vert (1))/w_res*2*pi;
a_boc = (vert(3) - vert (2))/w_res*2*pi;
a_cod = (vert(4) - vert (3))/w_res*2*pi;
a_doa = (vert(1) + w_res - vert (4))/w_res*2*pi;

% energy
v_ao  = [x0 y0]; v_bo = [x0-1 y0]; v_co = [x0-1 y0 - yc]; v_do = [x0 y0-yc];
n_v_ao = norm(v_ao); n_v_bo = norm(v_bo); n_v_co = norm(v_co); n_v_do = norm(v_do);
b_aob = acos(dot(v_ao, v_bo)/n_v_ao/n_v_bo);
if det([v_ao;v_bo]) < 0
    b_aob = 2*pi - b_aob;
end
b_boc = acos(dot(v_bo, v_co)/n_v_bo/n_v_co);
if det([v_bo;v_co]) < 0
    b_boc = 2*pi - b_boc;
end
b_cod = acos(dot(v_co, v_do)/n_v_co/n_v_do);
if det([v_co;v_do]) < 0
    b_cod = 2*pi - b_cod;
end
b_doa = acos(dot(v_do, v_ao)/n_v_do/n_v_ao);
if det([v_do;v_ao]) < 0
    b_doa = 2*pi - b_doa;
end

f = (b_aob - a_aob)^2 + (b_boc - a_boc)^2 + (b_cod - a_cod)^2 + (b_doa - a_doa)^2;

% gradient
% x0
d_aob_x0 = (2*x0-1)*n_v_ao*n_v_bo + dot(v_ao, v_bo) * (x0*n_v_bo/n_v_ao + (x0-1)*n_v_ao/n_v_bo);
d_aob_x0 = d_aob_x0 * (-1/(sqrt(1-cos(b_aob)*cos(b_aob))+eps))/n_v_ao/n_v_ao/n_v_bo/n_v_bo;
if det([v_ao;v_bo]) < 0
    d_aob_x0 = -d_aob_x0;
end
d_aob_x0 = 2*(b_aob - a_aob) * d_aob_x0;

d_boc_x0 = 2*(x0-1)*n_v_bo*n_v_co + dot(v_bo, v_co) * ((x0-1)*n_v_co/n_v_bo + (x0-1)*n_v_bo/n_v_co);
d_boc_x0 = d_boc_x0 * (-1/(sqrt(1-cos(b_boc)*cos(b_boc))+eps))/n_v_bo/n_v_bo/n_v_co/n_v_co;
if det([v_bo;v_co]) < 0
    d_boc_x0 = -d_boc_x0;
end
d_boc_x0 = 2*(b_boc - a_boc) * d_boc_x0;

d_cod_x0 = (2*x0-1)*n_v_co*n_v_do + dot(v_co, v_do) * ((x0-1)*n_v_do/n_v_co + x0*n_v_co/n_v_do);
d_cod_x0 = d_cod_x0 * (-1/(sqrt(1-cos(b_cod)*cos(b_cod))+eps))/n_v_co/n_v_co/n_v_do/n_v_do;
if det([v_co;v_do]) < 0
    d_cod_x0 = -d_cod_x0;
end
d_cod_x0 = 2*(b_cod - a_cod) * d_cod_x0;

d_doa_x0 = 2*x0*n_v_do*n_v_ao +dot(v_do, v_ao) * (x0*n_v_ao/n_v_do + x0*n_v_do/n_v_ao);
d_doa_x0 = d_doa_x0 * (-1/(sqrt(1-cos(b_doa)*cos(b_doa))+eps))/n_v_do/n_v_do/n_v_ao/n_v_ao;
if det([v_do;v_ao]) < 0
    d_doa_x0 = -d_doa_x0;
end
d_doa_x0 = 2*(b_doa - a_doa) * d_doa_x0;

% y0
d_aob_y0 = 2*y0*n_v_ao*n_v_bo + dot(v_ao, v_bo) * (y0*n_v_bo/n_v_ao + y0*n_v_ao/n_v_bo);
d_aob_y0 = d_aob_y0 * (-1/(sqrt(1-cos(b_aob)*cos(b_aob))+eps))/n_v_ao/n_v_ao/n_v_bo/n_v_bo;
if det([v_ao;v_bo]) < 0
    d_aob_y0 = -d_aob_y0;
end
d_aob_y0 = 2*(b_aob - a_aob) * d_aob_y0;

d_boc_y0 = (2*y0-yc)*n_v_bo*n_v_co + dot(v_bo, v_co) * (y0*n_v_co/n_v_bo + (y0-yc)*n_v_bo/n_v_co );
d_boc_y0 = d_boc_y0/n_v_bo/n_v_bo/n_v_co/n_v_co;
if det([v_bo;v_co]) < 0
    d_boc_y0 = -d_boc_y0;
end
d_boc_y0 = 2*(b_boc - a_boc) * d_boc_y0;

d_cod_y0 = 2*(y0-yc)*n_v_co*n_v_do + dot(v_co, v_do) * ((y0-yc)*n_v_do/n_v_co + (y0-yc)*n_v_co/n_v_do);
d_cod_y0 = d_cod_y0 * (-1/(sqrt(1-cos(b_cod)*cos(b_cod))+eps))/n_v_co/n_v_co/n_v_do/n_v_do;
if det([v_co;v_do]) < 0
    d_cod_y0 = -d_cod_y0;
end
d_cod_y0 = 2*(b_cod - a_cod) * d_cod_y0;

d_doa_y0 = (2*y0-yc)*n_v_do*n_v_ao + dot(v_do, v_ao) * ((y0-yc)*n_v_ao/n_v_do + y0*n_v_do/n_v_ao);
d_doa_y0 = d_doa_y0 * (-1/(sqrt(1-cos(b_doa)*cos(b_doa))+eps))/n_v_do/n_v_do/n_v_ao/n_v_ao;
if det([v_do;v_ao]) < 0
    d_doa_y0 = -d_doa_y0;
end
d_doa_y0 = 2*(b_doa - a_doa) * d_doa_y0;

% yc
d_boc_yc = (-y0)*n_v_bo*n_v_co + dot(v_bo, v_co) * n_v_bo* (yc-y0)/n_v_co;
d_boc_yc = d_boc_yc/n_v_bo/n_v_bo/n_v_co/n_v_co;
if det([v_bo;v_co]) < 0
    d_boc_yc = -d_boc_yc;
end
d_boc_yc = 2*(b_boc - a_boc) * d_boc_yc;

d_cod_yc = 2*(yc-y0)*n_v_co*n_v_do + dot(v_co, v_do) * ((yc-y0)*n_v_do/n_v_co+(yc-y0)*n_v_co/n_v_do);
d_cod_yc = d_cod_yc * (-1/(sqrt(1-cos(b_cod)*cos(b_cod))+eps))/n_v_co/n_v_co/n_v_do/n_v_do;
if det([v_co;v_do]) < 0
    d_cod_yc = -d_cod_yc;
end
d_cod_yc = 2*(b_cod - a_cod) * d_cod_yc;

d_doa_yc = (-y0)*n_v_do*n_v_ao + dot(v_do, v_ao) * n_v_ao * (yc-y0)/n_v_do;
d_doa_yc = d_doa_yc * (-1/(sqrt(1-cos(b_doa)*cos(b_doa))+eps))/n_v_do/n_v_do/n_v_ao/n_v_ao;
if det([v_do;v_ao]) < 0
    d_doa_yc = -d_doa_yc;
end
d_doa_yc = 2*(b_doa - a_doa) * d_doa_yc;

d_x0 = d_aob_x0 + d_boc_x0 + d_cod_x0 + d_doa_x0;
d_y0 = d_aob_y0 + d_boc_y0 + d_cod_y0 + d_doa_y0;
d_yc = d_boc_yc + d_cod_yc + d_doa_yc;

df = [d_x0; d_y0; d_yc; 0; 0; 0; 0; 0];

