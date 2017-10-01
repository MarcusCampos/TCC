
[baseName, folder] = uigetfile();
selectedFile = fullfile(folder, baseName);

I = imread(selectedFile);
figure = I
hold on, 
imshow(I);
[x,y,frame] = size (I);
whos I;
c = [1 12 146 410];
r = [1 104 156 129];
global x1;
global y1;
global x2;
global y2;
v = 0;
for j =1 :899
   if v == 0
    for i=1:1199
        redValue = I(j, i, 1);
        greenValue = I(j, i, 2);
        blueValue = I(j, i, 3);
        if redValue < 255 && greenValue <255 && blueValue == 255
            disp(redValue);
            disp(greenValue);
            disp(blueValue);
            x1 = i;
            y1 = j;
            v=1;
            break;
        end
    end 
    v=0;
     for i=1199:-1:1
            redValue2 = I(j, i, 1);
            greenValue2 = I(j, i, 2);
            blueValue2 = I(j, i, 3);
            if redValue2 == 0 && greenValue2 == 0 && blueValue2 == 255
                disp(redValue2);
                disp(greenValue2);
                disp(blueValue2);
                 x2 = i;
                y2 = j;
                v=1;
                break;
            end
     end    
   end
end 
v = 0;

plot(x1, y1);
plot(x2, y2);

hold off;
p= impixel(I,c,r)
