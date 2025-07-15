function fobfd_out=fobfd(x) 
% FOBFD is a function for fourth order backward finite difference
% this function was created because when doing numerical integration (ie
% using ODE23), the solver uses its in-built step sizes, which may not be
% constant. Hence this fourth order backward finite difference function was
% created to account for this. To fully understand it, some derivation is
% required.

a0=x(1);
a_1=x(2);
a_2=x(3);
a_3=x(4);
a_4=x(5);
f_1=x(7);
f_2=x(8);
f_3=x(9);
f_4=x(10);
df1n1=(f_1-f_2)/(a_1-a_2);
df1n2=(f_2-f_3)/(a_2-a_3);
df1n3=(f_3-f_4)/(a_3-a_4); 
df2n1=(df1n1-df1n2)/(a_1-a_2);
df2n2=(df1n2-df1n3)/(a_2-a_3);
df3n1=(df2n1-df2n2)/(a_1-a_2); 
df1n0=df1n1+df2n1*(a0-a_1)/factorial(2)+df3n1*(a0-a_1)^2/factorial(3);  
fobfd_out=df1n0;
end