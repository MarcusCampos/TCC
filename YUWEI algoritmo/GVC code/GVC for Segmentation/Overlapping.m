function N = Overlapping(figure,figure2)
I = imread(figure);
N = 0;
I2 = imread(figure2);
[x,y,frame] = size (I2);
x1 =0;
y1 =0;
x2 = 0;
y2 = 0;
intersec = 0;
uniao = 0; 
for i=1:y-1
    for j =1 :x-1
        redValue = I(j, i, 1);
        greenValue = I(j, i, 2);
        blueValue = I(j, i, 3);
        if (redValue < 100 && greenValue < 100 && blueValue < 100)
             redValue2 = I2(j, i, 1);
            greenValue2 = I2(j, i, 2);
            blueValue2 = I2(j, i, 3);
            if  (redValue == redValue2 && greenValue2== greenValue && blueValue == blueValue2)
                intersec = intersec +1;
                uniao = uniao +1;
            else 
                uniao=uniao+1;
            end
        end
    end 
end 
SepraStringPasta = strsplit(figure, '\');
ArquivoSaida = strsplit(figure, SepraStringPasta{14});
NomeArquivoResultados = strcat(ArquivoSaida{1}, 'ResultadosMetricos.txt');
fid = fopen(NomeArquivoResultados, 'at' );
overlapping = intersec/uniao;
%http://www.isprs.org/proceedings/XXXVII/congress/4_pdf/208.pdf
%https://clouard.users.greyc.fr/Pantheon/experiments/evaluation/index-en.html#example
undersegmentation = 1- intersec/ uniao;
oversegmentation = uniao+intersec / uniao;
FPR = uniao-intersec;
%combination measure ?
% OverSegmentationij = 1 - area(xi ? yj) / area(xi).
% UnderSegmentationij = 1 - area(xi ? yj) / area(yj). 
%https://www.researchgate.net/publication/280670971_Over-_and_Under-Segmentation_Evaluation_based_on_the_Segmentation_Covering_Measure
underlapping = intersec/intersec+uniao;
fprintf( fid, ' OverLapping do %s : %f\n', SepraStringPasta{15}, overlapping);
fclose(fid);
end
