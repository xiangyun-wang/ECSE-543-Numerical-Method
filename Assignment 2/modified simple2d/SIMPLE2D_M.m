function [Potential,Sglob] = SIMPLE2D_M(filename1,filename2,filename3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is solely developed to solve the ECSE 543 assignments.
% Send your questions and feedbacks to ali.akbarzadehsharbaf@mail.mcgill.ca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1
    % Open the input file
    fid = fopen(filename1,'r');
    data = textscan(fid, '%f%f%f%f',  'CollectOutput', 1);
    fclose(fid);
    data = data{1};
    
    % Extract the boundary condition data
    BC_nod = data(isnan(data(:,3)),1:2);
    
    % Extract x and y of each node
    Nod2xy = find(isnan(data(:,4)));
    Nod2xy = Nod2xy<(size(data,1)-size(BC_nod,1));
    Nod2xy = data(Nod2xy,1:3);
    Nod2xy = sortrows(Nod2xy,1);
    Nod2xy = transpose(Nod2xy(:,2:3));
    
    % Extract the nodes of each element
    Tri2Nod = transpose(data(  (size(Nod2xy,2)+1):((size(data,1)-size(BC_nod,1)))  ,1:3));
    size(Tri2Nod)
    Tri2Nod(:,1)
    Nod2xy(:,1)
    Nod2xy(:,2)
    Nod2xy(:,7)
elseif nargin == 3
    Nod2xy = dlmread(filename1);
    Nod2xy = transpose(Nod2xy(:,2:3));
    Tri2Nod = dlmread(filename2);
    Tri2Nod = Tri2Nod(:,1:3)';
    BC_nod = dlmread(filename3);
else
    error('Inappropriate Number of Inpt Files!')
end

% x and y corresponding to each node of every triangle
x1 = Nod2xy( 1 , Tri2Nod(1,:) );
x2 = Nod2xy( 1 , Tri2Nod(2,:) );
x3 = Nod2xy( 1 , Tri2Nod(3,:) );
y1 = Nod2xy( 2 , Tri2Nod(1,:) );
y2 = Nod2xy( 2 , Tri2Nod(2,:) );
y3 = Nod2xy( 2 , Tri2Nod(3,:) );

% Local stiffness matrix
b1 = y2 - y3;
b2 = y3 - y1;
b3 = y1 - y2;

c1 = x3 - x2;
c2 = x1 - x3;
c3 = x2 - x1;

Area = 0.5*abs(b1.*c2-b2.*c1);

Sloc = (0.25./repmat(Area,9,1)).*[ ...
    b1.*b1+c1.*c1 ; ...
    b1.*b2+c1.*c2 ; ...
    b1.*b3+c1.*c3 ; ...
    b1.*b2+c1.*c2 ; ...
    b2.*b2+c2.*c2 ; ...
    b3.*b2+c3.*c2 ; ...
    b1.*b3+c1.*c3 ; ...
    b2.*b3+c2.*c3 ; ...
    b3.*b3+c3.*c3 ];

% Assmble the local matrices into the global one
rows = Tri2Nod([1 2 3 1 2 3 1 2 3],:);
cols = Tri2Nod([1 1 1 2 2 2 3 3 3],:);
Sglob = sparse(rows,cols,Sloc);
b = zeros(length(Sglob),1); % the RHS vector

% Impose the Dirichlet BC (Boundary Condition)
Potential = SolveDirichlet(Sglob,b,BC_nod(:,1),BC_nod(:,2));

size(Potential)
size(Nod2xy)
Potential = [[1:length(Potential)]',Nod2xy',Potential];
pot = Potential(:,4);

% num = 0;
% 
% while(num < 46)
%     C(num*3+1,Tri2Nod(1,num+1)) = 1;
%     C(num*3+2,Tri2Nod(2,num+1)) = 1;
%     C(num*3+3,Tri2Nod(3,num+1)) = 1;
%     num = num + 1;
% end
% 
% num = 0;
% A_coefficient = 0.25 / (0.0002);
% while(num < 46)
%     x1 = Nod2xy(1,Tri2Nod(1,num+1));
%     x2 = Nod2xy(1,Tri2Nod(2,num+1));
%     x3 = Nod2xy(1,Tri2Nod(3,num+1));
%     y1 = Nod2xy(2,Tri2Nod(1,num+1));
%     y2 = Nod2xy(2,Tri2Nod(2,num+1));
%     y3 = Nod2xy(2,Tri2Nod(3,num+1));
%     Sdis(num*3 + 1, num*3 + 1) = A_coefficient * ((y2-y3)*(y2-y3)+(x3-x2)*(x3-x2));
%     Sdis(num*3 + 1, num*3 + 2) = A_coefficient * ((y2-y1)*(y3-y1)+(x3-x2)*(x1-x3));
%     Sdis(num*3 + 1, num*3 + 3) = A_coefficient * ((y2-y3)*(y1-y2)+(x3-x2)*(x2-x1));
%     Sdis(num*3 + 2, num*3 + 1) = A_coefficient * ((y2-y1)*(y3-y1)+(x3-x2)*(x1-x3));
%     Sdis(num*3 + 2, num*3 + 2) = A_coefficient * ((y3-y1)*(y3-y1)+(x1-x3)*(x1-x3));
%     Sdis(num*3 + 2, num*3 + 3) = A_coefficient * ((y3-y1)*(y1-y2)+(x1-x3)*(x2-x1));
%     Sdis(num*3 + 3, num*3 + 1) = A_coefficient * ((y2-y3)*(y1-y2)+(x3-x2)*(x2-x1));
%     Sdis(num*3 + 3, num*3 + 2) = A_coefficient * ((y3-y1)*(y1-y2)+(x1-x3)*(x2-x1));
%     Sdis(num*3 + 3, num*3 + 3) = A_coefficient * ((y1-y2)*(y1-y2)+(x2-x1)*(x2-x1));
% end

return


% Solve Ax=b subject to a given Dirichlet boundary condition
function X = SolveDirichlet(A,b,ConsVar,VarVal)

temp = sortrows([ConsVar VarVal],1);
ConsVar = temp(:,1);
VarVal = temp(:,2).';
Dim = length(b);
Free_var = setdiff(1:Dim,ConsVar);

A = A(Free_var,:);
b = b(Free_var,:);

b_mod = A(:,ConsVar).*repmat(VarVal,length(Free_var),1);
b_mod = sum(b_mod,2);
b = b - b_mod;
A = A(:,Free_var);

X_temp = A\b;

X=zeros(Dim,1);
X(ConsVar) = VarVal.';
X(Free_var) = X_temp;

return
