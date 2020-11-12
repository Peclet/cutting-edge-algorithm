function Proceed = PROCEED(Max_height,Min_height,Left_height, Right_height)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%This function decides if the data is worth analysing by checking the
%positions of the turning points on the profile.
%It terminates the analysis of a profile that is not worth pursuing
Benchmark = (0.6)*(Max_height);
Proceed = 1;

if Right_height < Benchmark
    Proceed=0;
elseif Left_height < Benchmark
    Proceed=0;
else 
end

end

