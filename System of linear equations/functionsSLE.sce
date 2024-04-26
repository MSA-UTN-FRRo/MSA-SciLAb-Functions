//Department of Chemical Engineering, Universidad Tecnológica Nacional - Facultad Regional Rosario, Argentina
//Course: Matemática Superior Aplicada
//Scilab functions for System of linear equations

function [D,e]=simple_gaussian(A,b)
    //Simple Gaussian elimination v.2024
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

function [C, P]=pivoteototal(C, i, P)
    // Universidad Tecnológica Nacional-Facultad Regional Rosario
    // Departamento de Ingeniería Química - Matemática Superior Aplicada
    // Alumno: Mignacco Mateo Leg: 51736  mateomignacco32@gmail.com
    //-------------------------------------------------------------
    //Función: Pivoteo total v.2024
    //-------------------------------------------------------------    
    [f c]=size(C);
    [big p]=max(abs(C([i:f],[i:f])));
    fp=p(1);
    cp=p(2);
    if fp > 1
        dummy = C;
        C(i,:) = dummy(i+fp-1,:);
        C(i+fp-1,:)= dummy(i,:);
    end
    if cp > 1
        dummy = C;
        C(:,i) = dummy(:,i+cp-1);
        C(:,i+cp-1)= dummy(:,i);
        dummy = P;
        P(:,i) = dummy(:,i+cp-1);
        P(:,i+cp-1)= dummy(:,i);
    end
endfunction

function C=pivoteoparcial(C, i)
    // Universidad Tecnológica Nacional-Facultad Regional Rosario
    // Departamento de Ingeniería Química - Matemática Superior Aplicada
    // Alumno: Mignacco Mateo Leg: 51736  mateomignacco32@gmail.com
    //-------------------------------------------------------------
    //Función: Pivoteo parcial v.2024
    //-------------------------------------------------------------    
    [f c]=size(C);
    [big p] = max(abs(C([i:f],i)));
    if p > 1
        dummy = C;
        C(i,:) = dummy(i+p-1,:);
        C(i+p-1,:)= dummy(i,:);
    end
endfunction

function [D, e]=gaussianaPP(A, b)
    // Universidad Tecnológica Nacional-Facultad Regional Rosario
    // Departamento de Ingeniería Química - Matemática Superior Aplicada
    // Alumno: Mignacco Mateo Leg: 51736  mateomignacco32@gmail.com
    //-------------------------------------------------------------
    //Función: Eliminacion gaussiana con pivoteo parcial v.2024
    //-------------------------------------------------------------
    [f c]=size(A);
    C=[A,b];
    for i=1:f-1
        C = pivoteoparcial(C, i)
        for j=i+1:f
            C(j,[i+1:c+1]) = C(j,[i+1:c+1])-(C(j,i)/C(i,i))*C(i,[i+1:c+1]);
        end
        C([i+1:f],i)=0;
    end
    D=C(:,[1:c]);
    e=C(:,c+1);
endfunction

function [D, e, P]=gaussianaPT(A, b)
    // Universidad Tecnológica Nacional-Facultad Regional Rosario
    // Departamento de Ingeniería Química - Matemática Superior Aplicada
    // Alumno: Mignacco Mateo Leg: 51736  mateomignacco32@gmail.com
    //-------------------------------------------------------------
    //Función: Eliminacion gaussiana con pivoteo total v.2024
    //-------------------------------------------------------------
    [f c]=size(A);
    C=[A,b];
    P = eye(f,f)
    for i=1:f-1
        [C,P] = pivoteototal(C, i, P)
        for j=i+1:f
            C(j,[i+1:c+1]) = C(j,[i+1:c+1])-(C(j,i)/C(i,i))*C(i,[i+1:c+1]);
        end
        C([i+1:f],i)=0;
    end
    D=C(:,[1:c]);
    e=C(:,c+1);
endfunction 
