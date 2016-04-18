function [seq]=MLseq(connec);

i=1;
j=1;
m=length(connec);
L=2^m -1;
registers=[zeros(1, m-1) 1];           %initial register contents
seq(i)=registers(m);                   %first element of the sequence
for i=2:L;
    new_reg_cont(1)=connec(1)*seq(i-1);
    for j=2:m,
        new_reg_cont(j)=registers(j-1) + connec(j)*seq(i-1);
    end;
    
    registers=new_reg_cont;            %current register content
    seq(i)=registers(m);               %next element of the sequence
end