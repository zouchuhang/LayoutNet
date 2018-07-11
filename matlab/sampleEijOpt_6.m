% 

function [f, df] = sampleEijOpt_6(x)
% optimization code for six wall case

% param
vert = x(6:11);
w_res = x(12);

x0 = x(1); y0 = x(2); % camera center
yc = x(3); % corner
xd = x(4);
yf = x(5);

% assume vertical lines are fixed
a_aob = (vert(2) - vert (1))/w_res*2*pi;
a_boc = (vert(3) - vert (2))/w_res*2*pi;
a_cod = (vert(4) - vert (3))/w_res*2*pi;
a_doe = (vert(5) - vert (4))/w_res*2*pi;
a_eof = (vert(6) - vert (5))/w_res*2*pi;
a_foa = (vert(1) + w_res - vert (6))/w_res*2*pi;

% energy
v_ao  = [x0 y0]; v_bo = [x0-1 y0]; v_co = [x0-1 y0 - yc]; 
v_do = [x0-xd y0-yc];  v_eo = [x0-xd y0-yf]; v_fo = [x0 y0-yf];

n_v_ao = norm(v_ao); n_v_bo = norm(v_bo); n_v_co = norm(v_co); 
n_v_do = norm(v_do); n_v_eo = norm(v_eo); n_v_fo = norm(v_fo);

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
b_doe = acos(dot(v_do, v_eo)/n_v_do/n_v_eo);
if det([v_do;v_eo]) < 0
    b_doe = 2*pi - b_doe;
end
b_eof = acos(dot(v_eo, v_fo)/n_v_eo/n_v_fo);
if det([v_eo;v_fo]) < 0
    b_eof = 2*pi - b_eof;
end
b_foa = acos(dot(v_fo, v_ao)/n_v_fo/n_v_ao);
if det([v_fo;v_ao]) < 0
    b_foa = 2*pi - b_foa;
end

f = (b_aob - a_aob)^2 + (b_boc - a_boc)^2 + (b_cod - a_cod)^2 + ...
    (b_doe - a_doe)^2 + (b_eof - a_eof)^2 + (b_foa - a_foa)^2;

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

d_cod_x0 = (2*x0-xd-1)*n_v_co*n_v_do + dot(v_co, v_do) * ((x0-1)*n_v_do/n_v_co + (x0-xd)*n_v_co/n_v_do);
d_cod_x0 = d_cod_x0 * (-1/(sqrt(1-cos(b_cod)*cos(b_cod))+eps))/n_v_co/n_v_co/n_v_do/n_v_do;
if det([v_co;v_do]) < 0
    d_cod_x0 = -d_cod_x0;
end
d_cod_x0 = 2*(b_cod - a_cod) * d_cod_x0;

d_doe_x0 = 2*(x0-xd)*n_v_do*n_v_eo + dot(v_do, v_eo) * ((x0-xd)*n_v_eo/n_v_do + (x0-xd)*n_v_do/n_v_eo);
d_doe_x0 = d_doe_x0 * (-1/(sqrt(1-cos(b_doe)*cos(b_doe))+eps))/n_v_do/n_v_do/n_v_eo/n_v_eo;
if det([v_do;v_eo]) < 0
    d_doe_x0 = -d_doe_x0;
end
d_doe_x0 = 2*(b_doe - a_doe) * d_doe_x0;

d_eof_x0 = (2*x0-xd)*n_v_eo*n_v_fo + dot(v_eo, v_fo) * ((x0-xd)*n_v_fo/n_v_eo + x0*n_v_eo/n_v_fo);
d_eof_x0 = d_eof_x0 * (-1/(sqrt(1-cos(b_eof)*cos(b_eof))+eps))/n_v_eo/n_v_eo/n_v_fo/n_v_fo;
if det([v_eo;v_fo]) < 0
    d_eof_x0 = -d_eof_x0;
end
d_eof_x0 = 2*(b_eof - a_eof) * d_eof_x0;

d_foa_x0 = 2*x0*n_v_fo*n_v_ao + dot(v_fo, v_ao) * (x0*n_v_ao/n_v_fo + x0*n_v_fo/n_v_ao);
d_foa_x0 = d_foa_x0 * (-1/(sqrt(1-cos(b_foa)*cos(b_foa))+eps))/n_v_fo/n_v_fo/n_v_ao/n_v_ao;
if det([v_fo;v_ao]) < 0
    d_foa_x0 = -d_foa_x0;
end
d_foa_x0 = 2*(b_foa - a_foa) * d_foa_x0;

% y0
d_aob_y0 = 2*y0*n_v_ao*n_v_bo + dot(v_ao, v_bo) * (y0*n_v_bo/n_v_ao + y0*n_v_ao/n_v_bo);
d_aob_y0 = d_aob_y0 * (-1/(sqrt(1-cos(b_aob)*cos(b_aob))+eps))/n_v_ao/n_v_ao/n_v_bo/n_v_bo;
if det([v_ao;v_bo]) < 0
    d_aob_y0 = -d_aob_y0;
end
d_aob_y0 = 2*(b_aob - a_aob) * d_aob_y0;

d_boc_y0 = (2*y0-yc)*n_v_bo*n_v_co + dot(v_bo, v_co) * (y0*n_v_co/n_v_bo + (y0-yc)*n_v_bo/n_v_co );
d_boc_y0 = d_boc_y0 * (-1/(sqrt(1-cos(b_boc)*cos(b_boc))+eps))/n_v_bo/n_v_bo/n_v_co/n_v_co;
if det([v_bo;v_co]) < 0
    d_boc_y0 = -d_boc_y0;
end
d_boc_y0 = 2*(b_boc - a_boc) * d_boc_y0;

d_cod_y0 = 2*(y0-yc)*n_v_co*n_v_do + dot(v_co, v_do) * ((y0-yc)/n_v_co*n_v_do + (y0-yc)/n_v_do*n_v_co);
d_cod_y0 = d_cod_y0 * (-1/(sqrt(1-cos(b_cod)*cos(b_cod))+eps))/n_v_co/n_v_co/n_v_do/n_v_do;
if det([v_co;v_do]) < 0
    d_cod_y0 = -d_cod_y0;
end
d_cod_y0 = 2*(b_cod - a_cod) * d_cod_y0;

d_doe_y0 = (2*y0-yf-yc)*n_v_do*n_v_eo + dot(v_do, v_eo)*((y0-yc)/n_v_do*n_v_eo + (y0-yf)/n_v_eo*n_v_do);
d_doe_y0 = d_doe_y0 * (-1/(sqrt(1-cos(b_doe)*cos(b_doe))+eps))/n_v_do/n_v_do/n_v_eo/n_v_eo;
if det([v_do;v_eo]) < 0
    d_doe_y0 = -d_doe_y0;
end
d_doe_y0 = 2*(b_doe - a_doe) * d_doe_y0;

d_eof_y0 = 2*(y0-yf)*n_v_eo*n_v_fo + dot(v_eo, v_fo)*((y0-yf)/n_v_eo*n_v_fo + (y0-yf)/n_v_fo*n_v_eo);
d_eof_y0 = d_eof_y0 * (-1/(sqrt(1-cos(b_eof)*cos(b_eof))+eps))/n_v_eo/n_v_eo/n_v_fo/n_v_fo;
if det([v_eo;v_fo]) < 0
    d_eof_y0 = -d_eof_y0;
end
d_eof_y0 = 2*(b_eof - a_eof) * d_eof_y0;

d_foa_y0 = (2*y0-yf)*n_v_fo*n_v_ao + dot(v_fo, v_ao)*((y0-yf)/n_v_fo*n_v_ao + y0/n_v_ao*n_v_fo);
d_foa_y0 = d_foa_y0 * (-1/(sqrt(1-cos(b_foa)*cos(b_foa))+eps))/n_v_fo/n_v_fo/n_v_ao/n_v_ao;
if det([v_fo;v_ao]) < 0
    d_foa_y0 = -d_foa_y0;
end
d_foa_y0 = 2*(b_foa - a_foa) * d_foa_y0;

% yc
d_boc_yc = (-y0)*n_v_bo*n_v_co + dot(v_bo, v_co) * n_v_bo* (yc-y0)/n_v_co;
d_boc_yc = d_boc_yc/n_v_bo/n_v_bo/n_v_co/n_v_co;
if det([v_bo;v_co]) < 0
    d_boc_yc = -d_boc_yc;
end
d_boc_yc = 2*(b_boc - a_boc) * d_boc_yc;

d_cod_yc = 2*(yc-y0)*n_v_co*n_v_do + dot(v_co, v_do)*((yc-y0)/n_v_co*n_v_do + (yc-y0)/n_v_do*n_v_co);
d_cod_yc = d_cod_yc * (-1/(sqrt(1-cos(b_cod)*cos(b_cod))+eps))/n_v_co/n_v_co/n_v_do/n_v_do;
if det([v_co;v_do]) < 0
    d_cod_yc = -d_cod_yc;
end
d_cod_yc = 2*(b_cod - a_cod) * d_cod_yc;

d_doe_yc = (yf-y0)*n_v_do*n_v_eo + dot(v_do, v_eo)*(yc-y0)/n_v_do*n_v_eo;
d_doe_yc = d_doe_yc * (-1/(sqrt(1-cos(b_doe)*cos(b_doe))+eps))/n_v_do/n_v_do/n_v_eo/n_v_eo;
if det([v_do;v_eo]) < 0
    d_doe_yc = -d_doe_yc;
end
d_doe_yc = 2*(b_doe - a_doe) * d_doe_yc;

% xd
d_cod_xd = (1-x0)*n_v_co*n_v_do + dot(v_co, v_do)*(xd-x0)*n_v_co/n_v_do;
d_cod_xd = d_cod_xd * (-1/(sqrt(1-cos(b_cod)*cos(b_cod))+eps))/n_v_co/n_v_co/n_v_do/n_v_do;
if det([v_co;v_do]) < 0
    d_cod_xd = -d_cod_xd;
end
d_cod_xd = 2*(b_cod - a_cod) * d_cod_xd;

d_doe_xd = 2*(xd-x0)*n_v_do*n_v_eo + dot(v_do, v_eo)*((xd-x0)/n_v_do*n_v_eo + (xd-x0)/n_v_eo*n_v_do);
d_doe_xd = d_doe_xd * (-1/(sqrt(1-cos(b_doe)*cos(b_doe))+eps))/n_v_do/n_v_do/n_v_eo/n_v_eo;
if det([v_do;v_eo]) < 0
    d_doe_xd = -d_doe_xd;
end
d_doe_xd = 2*(b_doe - a_doe) * d_doe_xd;

d_eof_xd = (-x0)*n_v_eo*n_v_fo + dot(v_eo, v_fo)*(xd-x0)/n_v_eo*n_v_fo;
d_eof_xd = d_eof_xd * (-1/(sqrt(1-cos(b_eof)*cos(b_eof))+eps))/n_v_eo/n_v_eo/n_v_fo/n_v_fo;
if det([v_eo;v_fo]) < 0
    d_eof_xd = -d_eof_xd;
end
d_eof_xd = 2*(b_eof - a_eof) * d_eof_xd;

% yf
d_doe_yf = (yc-y0)*n_v_do*n_v_eo + dot(v_do, v_eo)*(yf-y0)/n_v_eo*n_v_do;
d_doe_yf = d_doe_yf * (-1/(sqrt(1-cos(b_doe)*cos(b_doe))+eps))/n_v_do/n_v_do/n_v_eo/n_v_eo;
if det([v_do;v_eo]) < 0
    d_doe_yf = -d_doe_yf;
end
d_doe_yf = 2*(b_doe - a_doe) * d_doe_yf;

d_eof_yf = 2*(yf-y0)*n_v_eo*n_v_fo + dot(v_eo, v_fo)*((yf-y0)/n_v_eo*n_v_fo + (yf-y0)/n_v_fo*n_v_eo);
d_eof_yf = d_eof_yf * (-1/(sqrt(1-cos(b_eof)*cos(b_eof))+eps))/n_v_eo/n_v_eo/n_v_fo/n_v_fo;
if det([v_eo;v_fo]) < 0
    d_eof_yf = -d_eof_yf;
end
d_eof_yf = 2*(b_eof - a_eof) * d_eof_yf;

d_foa_yf = (-y0)*n_v_fo*n_v_ao + dot(v_fo, v_ao)*(yf-y0)/n_v_fo*n_v_ao;
d_foa_yf = d_foa_yf * (-1/(sqrt(1-cos(b_foa)*cos(b_foa))+eps))/n_v_fo/n_v_fo/n_v_ao/n_v_ao;
if det([v_fo;v_ao]) < 0
    d_foa_yf = -d_foa_yf;
end
d_foa_yf = 2*(b_foa - a_foa) * d_foa_yf;

d_x0 = d_aob_x0 + d_boc_x0 + d_cod_x0 + d_doe_x0 + d_eof_x0 + d_foa_x0;
d_y0 = d_aob_y0 + d_boc_y0 + d_cod_y0 + d_doe_y0 + d_eof_y0 + d_foa_y0;
d_yc = d_boc_yc + d_cod_yc + d_doe_yc;
d_xd = d_cod_xd + d_doe_xd + d_eof_xd;
d_yf = d_doe_yf + d_eof_yf + d_foa_yf;

df = [d_x0; d_y0; d_yc; d_xd; d_yf; 0; 0; 0; 0; 0;0;0];

