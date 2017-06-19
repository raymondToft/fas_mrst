perm_range = [0.01 0.5]
G = cartGrid([48, 48, 4]);
filter_size = [3 3 3];
std_ = 2.5;
p= gaussianField(G.cartDims, perm_range, filter_size, std_);
        
K = p.^3. * (1e-5)^2./(0.81 * 72 * (1-p).^2);

maxK = max(max(max(K)))
minK = min(min(min(K)))
