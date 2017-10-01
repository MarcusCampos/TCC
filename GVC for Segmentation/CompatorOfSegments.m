function CompatorOfSegments(endox, endoy,epix, epiy,segendox, segendoy, segepix,segepiy)
sizeEndox = size(endox);
sizeEndoy  = size(endoy);
sizeEpix  = size(epix);
sizeEpiy  = size(epiy);
sizeSegendox  = size(segendox);
sizeSegendoy  = size(segendoy);
sizeSegepix = size(segepix);
sizeSegepiy = size(segepiy);

fileID = fopen('Results.txt','wt');
fprintf(fileID,'Arquivo para comprar os resultados obtidos');

if (sizeEndox(1,1) < sizeSegendox(1,1))
    fprintf(fileID,'Diferencia do número de pontos do endocCárdio em X,  %12.8f', sizeSegendox(1,1) - sizeEndox(1,1));
elseif  (sizeEndox(1,1) > sizeSegendox(1,1))
     fprintf(fileID,'Diferencia do número de pontos do endocCárdio em X,  %12.8f', sizeEndox(1,1) - sizeSegendox(1,1));
    else
     fprintf(fileID,'Não Há diferença de pontos');
end

if (sizeEndoy(1,1) < sizeSegendoy(1,1))
     fprintf(fileID,'Diferencia do número de pontos do endocCárdio em Y,  %12.8f', sizeSegendoy(1,1) - sizeEndoy(1,1));
 elseif  (sizeEndoy(1,1) > sizeSegendoy(1,1))
      fprintf(fileID,'Diferencia do número de pontos do endocCárdio em Y,  %12.8f', sizeEndoy(1,1) - sizeSegendoy(1,1));
    else
     fprintf(fileID,'Não Há diferença de pontos');
 end

if (sizeEpix(1,1) < sizeSegepix(1,1))
     fprintf(fileID,'Diferencia do número de pontos do epiCárdio em X,  %12.8f', sizeSegepix(1,1) - sizeEpix(1,1));
    elseif  (sizeEpix(1,1) > sizeSegepix(1,1))
         fprintf(fileID,'Diferencia do número de pontos do epiCárdio em X,  %12.8f', sizeEpix(1,1) - sizeSegepix(1,1));
    else
         fprintf(fileID,'Não Há diferença de pontos');
end

if (sizeEpiy(1,1) < sizeSegepiy(1,1))
     fprintf(fileID,'Diferencia do número de pontos do epiCárdio em Y,  %12.8f', sizeSegepiy(1,1) - sizeEpiy(1,1));
elseif  (sizeSegendox(1,1) > sizeSegepiy(1,1))
     fprintf(fileID,'Diferencia do número de pontos do epiCárdio em Y,  %12.8f', sizeEpiy(1,1) - sizeSegepiy(1,1));
    else
     fprintf(fileID,'Não Há diferença de pontos');
end

control = 0; 

for count = 1:sizeEndox
    for countseg = 1:sizeSegendox
        if endox(count,1) - segendox(countseg,1)< 1
            control = control +1;
        end
    end
end



fclose(fileID);

end