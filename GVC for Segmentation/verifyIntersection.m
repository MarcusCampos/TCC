function N = verifyIntersection(xI1, yI1, centerXI1, centerYI1, xI2, yI2, centerXI2, centerYI2)
     intersec = 0;
for i =1:size(xI1)
        for k =1 : size(yI1)
           if k == i
               A = xI1(i,1) ;
               B = yI1(k,1);
               for kl = A:max(xI1)
                   for o=1: size(xI1)
                       if A >= xI1(o,1)
                           break;
                       else 
                            for i2 =1:size(xI2)
                                 disp('Entrou '+ i2);
                            for k2 =1 : size(yI2)
                                  disp('Entrou '+ k2);
                               if k2 == i2
                                   A2 = xI2(i2,1) ;
                                   B2 = xI2(k2,1);
                                   for klB = A2:max(xI2)
                                       for o2=1: size(xI2)
                                           if A2 >= xI2(o2,1)
                                               break;
                                           else 
                                                 if round(A) == round(A2) && round(B)==round(B2);
                                                     intersec = intersec+1;
                                                     disp('Entrou '+ A + ' '+ B );
                                                 end
                                                 break;
                                           end
                                       end
                                       A2=A2+1;
                                   end
                               end 
                            end
                         end
                       end
                   end
                   A=A+1;
               end
           end 
        end
end
    N = intersec
     return 
end
