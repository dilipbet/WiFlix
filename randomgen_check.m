%clear all
%close all
    stream = RandStream('mt19937ar','Seed',28);
    RandStream.setDefaultStream(stream);
    for i = 1:1
        rand(1,1)
        poissrnd(1000,1,1)
    end
