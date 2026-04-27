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
// Alumno: Mignacco Mateo Leg: 51736
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
// Alumno: Mignacco Mateo Leg: 51736
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
// Alumno: Mignacco Mateo Leg: 51736
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

function x = SD(A,b)
    // Universidad Tecnológica Nacional-Facultad Regional Rosario
    // Departamento de Ingeniería Química - Matemática Superior Aplicada
    // Alumno: Liang Martina  Leg: 53728
    //-------------------------------------------------------------
    //Función: Sustitución hacia adelante v.2026
    //-------------------------------------------------------------
    [n c] = size(A)
    x(1,1) = b(1)/A(1,1)
    for i = 2:n
        x(i,1) = (b(i)-A(i,[1:i-1])*x([1:i-1],1))/A(i,i)
    end
endfunction

function x = SA(A,b)
    // Universidad Tecnológica Nacional-Facultad Regional Rosario
    // Departamento de Ingeniería Química - Matemática Superior Aplicada
    // Alumno: Liang Martina  Leg: 53728
    //-------------------------------------------------------------
    //Función: Sustitución hacia atrás v.2026
    //-------------------------------------------------------------
    [n c] = size(A)
    x(n,1) = b(n)/A(n,n)
    for i=n-1:-1:1
        x(i,1) = (b(i)-A(i,[i+1:n])*x([i+1:n],1))/A(i,i)
    end
endfunction

function AINV = MYINV(A)
    // Universidad Tecnológica Nacional-Facultad Regional Rosario
    // Departamento de Ingeniería Química - Matemática Superior Aplicada
    // Alumno: Liang Martina  Leg: 53728
    //-------------------------------------------------------------
    //Función: Inversa de una matriz v.2026
    //-------------------------------------------------------------
    [L U P] = lu(A)
    [n c] = size(A)
    AINV = []
    Id = eye(n,c)
    for i = 1:c
        b = Id(:,i)
        z = P * b
        y = SD(L,z)
        x = SA(U,y)
        AINV(:,i) = x
    end
endfunction


function x = Thomas(M,b)
    // Universidad Tecnológica Nacional-Facultad Regional Rosario
    // Departamento de Ingeniería Química - Matemática Superior Aplicada
    // Alumno: Liang Martina  Leg: 53728
    //-------------------------------------------------------------
    //Función: Thomas v.2026
    //-------------------------------------------------------------
    A = diag(M,-1)
    B = diag(M,0)
    C = diag(M,1)
    A = [0;A]
    C = [C;0]
    
    n = length(b)
    
    P = []
    Q = []
    x = []
    P(1) = C(1)/B(1)
    Q(1) = b(1)/B(1)
    
    for i=2:(n-1)
        P(i) = C(i)/(B(i)-A(i)*P(i-1))
        Q(i) = (b(i)-A(i)*Q(i-1))/(B(i)-A(i)*P(i-1))
    end
    
    Q(n) = (b(n)-A(n)*Q(n-1))/(B(n)-A(n)*P(n-1))
    x(n) = Q(n)
    
    for i=(n-1):-1:1
        x(i) = Q(i)-P(i)*x(i+1)
    end
endfunction
