function [T1_X,T1_Y,T2_X, T2_Y] = TurningPoints(XY1,Y,S2,n)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
counter=1;
for j=1
    for i = 1:((n-1))
        if S2(i+1,1) == S2(i,1)
        else
        idxr(counter,1)=i;
        counter = counter+1;
        end
    end
    if size(idxr,1)==3
        TurnLoc = zeros(3,2);
        y=3;
        TurnLoc(:,1) = XY1(idxr);
        TurnLoc(:,2) = Y(idxr);
        T1_X = TurnLoc(1,1);
        T1_Y = TurnLoc(1,2);
        T2_X = TurnLoc(3,1);
        T2_Y = TurnLoc(3,2);
    else
        [TF,S1,S2] = ischange(Y,'linear','MaxNumChanges',4);
        counter=1;     
        y=4;
        for i = 1:((n-1))
            if S2(i+1,1) == S2(i,1)
            else
                idxr(counter,1)=i;
                counter = counter+1;
            end
        end
        TurnLoc = zeros(4,2);
        TurnLoc(:,1) = XY1(idxr);
        TurnLoc(:,2) = Y(idxr);
        T1_X = TurnLoc(2,1);
        T1_Y = TurnLoc(2,2);
        T2_X = TurnLoc(4,1);
        T2_Y = TurnLoc(4,2);
    end 
end 

end

