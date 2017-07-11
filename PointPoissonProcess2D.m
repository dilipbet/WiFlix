function[Locs, NumPoints] =  PointPoissonProcess2D(cell_area,density, frame_variant, cell_side, NumFrames)


%[pdf,cdf] = poisson_pdf(density*cell_area);


%NumPoints= density*cell_area*ones(NumFrames,1); %

% NumPoints= poisson_sample(density*cell_area,NumFrames);
  NumPoints = poissrnd(density*cell_area,1,1);

if(frame_variant == 1)
    %each row is for a new frame
    %create MaxNumPoints many points in the square region
    Locs = (rand(NumFrames,max(NumPoints))+rand(NumFrames,max(NumPoints))*sqrt(-1))*cell_side;
elseif(frame_variant ==0)
    % create according to the first row and keep their locations and
    % numbers fixed
    Locs = (rand(1,max(NumPoints(1)))+rand(1,max(NumPoints(1)))*sqrt(-1))*cell_side;
    Locs = repmat(Locs,NumFrames,1);
    NumPoints = repmat(NumPoints(1),NumFrames);

end


