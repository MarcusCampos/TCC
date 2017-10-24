function [ P1,P2,P3 ] = getPontosIniciais( stringFilepath )
    I = imread(stringFilepath);
    N = 0;
    Fig = I;
    hold on
    [x,y,frame] = size (I);
    x1 =0;
    y1 =0;
    x2 = 0;
    y2 = 0;
    pontos = cell(3,2);
    fig = figure(1); 
    hold off;
    figure(1)
    imshow(I);
    hold on,
    for i=1:y-1

        for j =1 :x-1
            redValue = I(j, i, 1);
            greenValue = I(j, i, 2);
            blueValue = I(j, i, 3);
            if greenValue >= 240
                disp('oi');
            end
            if redValue < 100 && blueValue<100 &&greenValue >= 200
                disp('tchau');
            end
            if (redValue < 100 && blueValue<100 &&greenValue >= 200) 
                disp(redValue);
                disp(greenValue);
                disp(blueValue);
                x1 = i;
                y1 = j;
                 if isempty(pontos(1,1)) ==0 
                    pontos(1,1) = {i};
                    pontos(1,2) = {j};
                    break;
                end
                if isempty(pontos(2,1)) == 0
                    pontos(2,1) = {i};
                    pontos(2,2) = {j};
                    break;
                else
                      pontos(3,1) = {i};
                        pontos(3,2) = {j};
                        break;
                end;
           
            end
        end 
    end 
     pontos(1), pontos(2), pontos(3);
end

