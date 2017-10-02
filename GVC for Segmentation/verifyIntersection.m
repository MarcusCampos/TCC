function N = verifyIntersection(figure,figure2)
I = imread(figure);
N = 0;
I2 = imread(figure2);
Fig = I
hold on
imshow(I);
[x,y,frame] = size (I);
x1 =0;
y1 =0;
x2 = 0;
y2 = 0;

hold on
for i=1:1199

    for j =1 :898
        redValue = I(j, i, 1);
        greenValue = I(j, i, 2);
        blueValue = I(j, i, 3);
        if (redValue < 100 && greenValue < 100 && blueValue < 100) 
            disp(redValue);
            disp(greenValue);
            disp(blueValue);
            x1 = i;
            y1 = j;
           
            break;
        end
    end 
    for n=898:-1:1
            redValue2 = I(n, i, 1);
            greenValue2 = I(n, i, 2);
            blueValue2 = I(n, i, 3);
            if (redValue2 <100 && greenValue2 < 100 && blueValue2 < 100)  
                disp(redValue2);
                disp(greenValue2);
                disp(blueValue2);
                 x2 = i;
                y2 = n;
                plot([x2,x1],[y2,y1],'k','LineWidth',1);
                break;
           end
    end  
 
end 
hold off;
strgeral = 'Preenchido1';
fileNameFormat = '.jpg';
stringNameFile = strcat(strgeral,fileNameFormat);
saveas(Fig, figure+stringNameFile);

Fig2 = I2
imshow(I2);
[x,y,frame] = size (I2);
x1 =0;
y1 =0;
x2 = 0;
y2 = 0;

hold on
for i=1:1199

    for j =1 :898
        redValue = I2(j, i, 1);
        greenValue = I2(j, i, 2);
        blueValue = I2(j, i, 3);
        if (redValue < 100 && greenValue < 100 && blueValue < 100) 
            disp(redValue);
            disp(greenValue);
            disp(blueValue);
            x1 = i;
            y1 = j;
           
            break;
        end
    end 
    for n=898:-1:1
            redValue2 = I2(n, i, 1);
            greenValue2 = I2(n, i, 2);
            blueValue2 = I2(n, i, 3);
            if (redValue2 <100 && greenValue2 < 100 && blueValue2 < 100)  
                disp(redValue2);
                disp(greenValue2);
                disp(blueValue2);
                 x2 = i;
                y2 = n;
                plot([x2,x1],[y2,y1],'k','LineWidth',1);
                break;
           end
    end  
 
end 
hold off;
strgeral = 'Preenchido2';
stringnumer = num2str(k);
fileNameFormat = '.jpg';
stringNameFile = strcat(strgeral,fileNameFormat);
saveas(Fig2, figure2+stringNameFile);

I = imread(figure);
I2 = imread(figure2);

for i=1:1199
    for j =1 :898
        redValue = I(j, i, 1);
        greenValue = I(j, i, 2);
        blueValue = I(j, i, 3);
        if (redValue < 100 && greenValue < 100 && blueValue < 100)
             redValue2 = I2(j, i, 1);
            greenValue2 = I2(j, i, 2);
            blueValue2 = I2(j, i, 3);
            if  (redValue == redValue2 && greenValue2== greenValue && blueValue == blueValue2)
                N=N+1;
            end
        end
    end 
end 
end
