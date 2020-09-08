

%**********Specify the constants*************%
mu= 10;
kappa= 1.3;
wa= 0.3;
a= 45;

%********Calculate the Right Cauchy Green Tensor******%

syms CC11 CC12 CC21 CC22 real

C(1,1) = CC11;
C(1,2) = CC12;
C(2,1) = CC21;
C(2,2) = CC22;
C(3,3) = 1/(CC11*CC22 - CC12*CC21);
C(1,3) = 0;
C(3,1) = 0;
C(2,3) = 0;
C(3,2) = 0;

C_sym = 0.5*(C + transpose(C));

%********Defining the structural tensor and pressure term*********%

L = [(wa*cosd(a)^2)+(1-wa)/3 wa*cosd(a)*sind(a) 0; wa*cosd(a)*sind(a) (wa*sind(a)^2)+(1-wa)/3 0; 0 0 (1-wa)/3 ];

pres = ((kappa*2*(1-wa)/3)*(1/(C_sym(1,1)*C_sym(2,2)-C_sym(1,2)^2)) + (mu*2*(wa-1)/3)*((C_sym(1,1)*C_sym(2,2)-C_sym(1,2)^2)));

%**********Equation for strain energy function********%

psi = kappa*(trace(C_sym*L)-1) + mu*(trace(inv(C_sym)*L)-1)- pres*(det(C_sym)-1);

%**********Equation for Piola Kirchoff Stress**********%

for i=1:2
    for j=1:2
        S(i,j) = 2*diff(psi,C(i,j));
    end
end

%**********Assumption of values to cross check with fortran**********%

CC11 = 17.277320000000003;
CC12 = 27.643700000000003;
CC21 = 27.643700000000003;
CC22 = 44.373620000000003;

format long

   eval(S)





