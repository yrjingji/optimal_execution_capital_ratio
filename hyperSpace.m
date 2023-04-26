sigmaSeq = linspace(0.01, 0.91, 5);
% ... and so on for all hyperparameters
initialpriceSeq = linspace(16,20,5);
tauSeq = linspace(0.1,1,5);
lSeq = linspace(800,1200,5);
x0Seq = linspace(250,350,5);
s0Seq = linspace(200,250,5);
betaSeq = linspace(0.08,0.09,5);
gammaSeq = linspace(0.06,0.08,5);
vectors = {sigmaSeq,initialpriceSeq,tauSeq,lSeq,x0Seq,s0Seq,betaSeq,gammaSeq};
all_combinations = get_combinations(vectors);