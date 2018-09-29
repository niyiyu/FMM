load v_origin.txt
v=imresize(v_origin,0.5,'bilinear');
dlmwrite('/Users/niyiyu/documents/fortran/fmm/main/v.txt',v,' ');
clear all