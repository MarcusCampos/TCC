function [ P1,P2,P3 ] = getPontosIniciais( stringFilepath )
    I = imread(stringFilepath)
    N = 0;
%     [x,y,frame] = size (I);
%     if x > 460   || y > 460
%         I = imresize(I,[460 460])
%     end
%     clearvars x y;
%     [x,y,frame] = size (I);
%     if y > 1800
%         I = imcrop(I, [714.5 0.5 452 456]);
%     end;
    [x,y,frame] = size (I);
    x1 =0;
    y1 =0;
    x2 = 0;
    y2 = 0;
    in = 1
    pontos = cell(3,2);
   % imshow(I);
    i=1;
     Fig = I;
    in =1;
    while in< x-1 
        i =1;
        while i<y-1
            redValue = Fig(in, i, 1);
            greenValue = Fig(in, i, 2);
            blueValue = Fig(in, i, 3);
            if (redValue ==255 && blueValue < 10 && greenValue <10)
                disp(redValue);
                disp(greenValue);
                disp(blueValue);
                x1 = i;
                y1 = in;
                 if isempty(pontos{1,1}) 
                    pontos(1,1) = {i};
                    pontos(1,2) = {in};
                    break;
                elseif isempty(pontos{2,1}) 
                    pontos(2,1) = {i};
                    pontos(2,2) = {in};
                    break;
                else
                      pontos(3,1) = {i};
                        pontos(3,2) = {in};
                   break;
                end;
           
            end;
            i = i +1 ;
        end 
        in = in +1 ;
    end 
    P1 = [pontos(2,1), pontos(2,2)];
    P2 = [pontos(3,1),pontos(3,2)]  ;
   P3 = [pontos(1,1), pontos(1,2)];
end

