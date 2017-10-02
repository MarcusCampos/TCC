
[baseName, folder] = uigetfile();
selectedFile = fullfile(folder, baseName);

I = imread(selectedFile);
figure = I
hold on
imshow(figure);
[x,y,frame] = size (figure);
x1 =0;
y1 =0;
global x2;
global y2;

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
