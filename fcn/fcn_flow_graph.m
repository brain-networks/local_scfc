function dyn = fcn_flow_graph(aij,r,t)
%FCN_flow_graph       continuous time random walk on network
%
%   DYN = FCN_flow_graph(AIJ,R,T) instantiates a continuous time
%         random walk on network with adjacency matrix AIJ. Waiting time
%         for walkers at each node are distributed as Poisson with rate
%         parameter R. This function returns the flow graph at time T.
%
%   Inputs:         AIJ,        symmetric adjacency matrix
%                     R,        rate parameter
%                     T,        markov time
%
%   Outputs:        DYN,        flow graph at time T
%
%   Reference: Lambiotte et al 2011. DOI: 10.1103/PhysRevE.84.017102
%
%   Richard Betzel, Indiana University, 2012
%

%modification history
%04.18.2012 - original
%07.20.2012 - added symmetrizing step (dyn + dyn')/2 to correct for
%             precision issues

[rr,cr] = size(r);
if rr > cr
    r = r';
end
s = sum(aij);
z = sum(s./r);
ps = s./(z.*r);
lap = -(bsxfun(@times,bsxfun(@rdivide,aij,s),r) - diag(r));

dyn = bsxfun(@times,z*expm(-t*lap),ps);
dyn = (dyn + dyn')./2;