function funtion_plots_shadings(xvec,minB,maxB,fillcolor)        

if isrow(xvec)==0
    xvec=xvec';
end

if isrow(minB)==0
    minB=minB';
end

if isrow(maxB)==0
    maxB=maxB';
end

    x22 = [xvec, fliplr(xvec)];
    inBetween2 = [minB, fliplr(maxB)];
    fill(x22, inBetween2,fillcolor,'LineStyle','none');
    alpha(0.05)
        
end