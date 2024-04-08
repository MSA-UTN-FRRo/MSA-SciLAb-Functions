//Department of Chemical Engineering, Universidad Tecnologica Nacional - Facultad Regional Rosario
//Course: Matemática Superior Aplicada
//Function: Simple Gaussian elimination
function [D,e]=gaussiana(A,b)
    [f c]=size(A);
    C=[A,b];
    for i=1:f-1
        for j=i+1:f
            C(j,:)=C(j,:)-(C(j,i)/C(i,i))*C(i,:);
        end
    end
    D=C(:,[1:c]);
    e=C(:,c+1);
endfunction
