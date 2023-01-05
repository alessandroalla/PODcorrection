switch choice
    case 0 % Time interval [0,T]
        if model == 1
            R = 45;
            rdeim = 48;
        elseif model == 2
            R = 80;
            rdeim = 90;
        else
            R = 200;
            rdeim = 363;
        end
    case 1 % Zone I1 reactivity
        if model == 1
            R = 45;
            rdeim = 50;
        elseif model == 2
            R = 35;
            rdeim = 50;
        else
            R = 41;
            rdeim = 60;
        end
    case 2 % Zone I2
        if model == 2
            R = 15;
            rdeim = 14;
        elseif model == 3
            R = 150;
            rdeim = 324;
        else
            rdeim = 13;
            R = 12;
        end
end