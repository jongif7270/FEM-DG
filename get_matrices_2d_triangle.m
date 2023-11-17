function [M_R,Srr_R,Srs_R,Ssr_R,Sss_R,Dr_R,Ds_R]=get_matrices_2d_triangle(k)
if k==1
M_R = [2 1 1; 1 2 1; 1 1 2]/6;
Srr_R = [1 -1 0; -1 1 0; 0 0 0]/2;
Srs_R = [1 0 -1; -1 0 1; 0 0 0]/2;
Ssr_R = [1 -1 0; 0 0 0; -1 1 0]/2;
Sss_R = [1 0 -1; 0 0 0; -1 0 1]/2;
Dr_R = [-1 1 0; -1 1 0; -1 1 0]/2;
Ds_R = [-1 0 1; -1 0 1; -1 0 1]/2;
elseif k==2
M_R = [6 0 -1 0 -4 -1; 0 32 0 16 16 -4; -1 0 6 -4 0 -1;
0 16 -4 32 16 0; -4 16 0 16 32 0; -1 -4 -1 0 0 6]/90;
Srr_R = [3 -4 1 0 0 0; -4 8 -4 0 0 0; 1 -4 3 0 0 0;
0 0 0 8 -8 0; 0 0 0 -8 8 0; 0 0 0 0 0 0]/6;
Srs_R = [3 0 0 -4 0 1; -4 4 0 4 -4 0; 1 -4 0 0 4 -1;
0 4 0 4 -4 -4; 0 -4 0 -4 4 4; 0 0 0 0 0 0]/6;

Ssr_R = [3 -4 1 0 0 0; 0 4 -4 4 -4 0; 0 0 0 0 0 0;
-4 4 0 4 -4 0; 0 -4 4 -4 4 0; 1 0 -1 -4 4 0]/6;
Sss_R = [3 0 0 -4 0 1; 0 8 0 0 -8 0; 0 0 0 0 0 0;
-4 0 0 8 0 -4; 0 -8 0 0 8 0; 1 0 0 -4 0 3]/6;
Dr_R = [-3 4 -1 0 0 0; -1 0 1 0 0 0; 1 -4 3 0 0 0;
-1 2 -1 -2 2 0; 1 -2 1 -2 2 0; 1 0 -1 -4 4 0]/2;
Ds_R = [-3 0 0 4 0 -1; -1 -2 0 2 2 -1; 1 -4 0 0 4 -1;
-1 0 0 0 0 1; 1 -2 0 -2 2 1; 1 0 0 -4 0 3]/2;
else
M_R = 0; Srr_R = 0; Dr_R = 0;
end
end