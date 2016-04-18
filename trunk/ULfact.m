function [Upp] = ULfact(sq1,R)        % R=UL factorisation
for i=sq1:-1:1,
    for w=sq1:-1:1,
        temp=0;
        if (i==w),
            if (i<sq1),
                for n2=i+1:sq1,
                    temp = temp + Upp(i,n2)*Upp(i,n2);
                end;
            end;
            Upp(i,i)=sqrt(R(i,i)-temp);
        else
            if (w<sq1),
                for n2=w+1:sq1,
                    temp=temp + Upp(i,n2)*Upp(w,n2);
                end;
            end;
            if (w<i),
                Upp(i,w)=0;
            else
                Upp(i,w)=(R(i,w)- temp)/Upp(w,w);
            end;
        end;
    end;
end;
