//Department of Chemical Engineering, Universidad Tecnológica Nacional - Facultad Regional Rosario, Argentina
//Course: Matemática Superior Aplicada
//Scilab functions for System of linear equation

function [D,e]=simple_gaussian(A,b)
//Simple Gaussian elimination
    [f c]=size(A);
    C=[A,b];
    for i=1:f-1
        for j=i+1:f
            C(j,[i+1:c+1]) = C(j,[i+1:c+1])-(C(j,i)/C(i,i))*C(i,[i+1:c+1]);
        end
        C([i+1:f],i)=0;
    end
    D=C(:,[1:c]);
    e=C(:,c+1);
endfunction
