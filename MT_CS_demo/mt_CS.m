function [weights,ML] = mt_CS(PHI,t,a,b,eta)
%---------------------------------------------------------------------------------------
% The MT-CS algorithm for the following paper:
% "Multi-Task Compressive Sesning" (Preprint, 2007). The algorithm 
% is an extension of the fast RVM algorithm [Tipping & Faul, 2003]
% in two-fold: (i) the noise variance is marginalized and (ii) it is for
% multi-task CS, including single-task CS as a special case
% Coded by: Shihao Ji, ECE, Duke University
% last change: May. 15, 2007
% A bug was fixed on Aug.03, 2008 for the cases where signals are dramatic undersampled
%---------------------------------------------------------------------------------------
% Input:
%   PHI: projection matrix. Cell structure, One cell for one task.
%   t:   CS measurements. Cell structure, One cell for one task.
%   a,b: parameters of Gamma prior on noise variance
%   eta: threshold for stopping the algorithm (suggested value: 1e-8)
% Output:
%   weights: sparse weights for all the tasks. One column for one task
%   ML:      the increase of the joint mariginal likelihood for each
%            iteration
%
if iscell(t)
    NT = length(t);
else
    NT = 1;
    PHI = {PHI};
    t = {t};
end
fprintf(1,'This is a %d-task learning!\n',NT);
%
% find initial alpha
for k = 1:NT
    [N(k),M(k)] = size(PHI{k});
end
if sum(abs(M-M(1))) ~= 0
    error('Sorry! The sizes of the underlying signals should be the same!\n');
else
    M = M(1);
end
% find initial alpha
K = repmat(N+2*a,[M,1]);
%
for k = 1:NT
    PHIt(:,k) = PHI{k}'*t{k};
    PHI2(:,k) = sum(PHI{k}.^2)';
    G2(k) = t{k}'*t{k}+2*b;
end
G2 = repmat(G2,[M,1]);
X = G2.*PHI2./PHIt.^2;
ml = K.*log(X./K)-(K-1).*log((X-1)./(K-1));

ml_sum = sum(ml,2);
while 1
    [ML,index] = max(ml_sum);
    alpha = NT./sum((K(index,:).*PHIt(index,:).^2./G2(index,:)-PHI2(index,:))./(PHI2(index,:).*(PHI2(index,:)-PHIt(index,:).^2./G2(index,:))),2);
    if alpha > 0  % A bug was fixed here, Aug 03, 2008 (alpha should be greater than 0)
        break;
    else
        ml_sum(index) = 0;
    end
end

for k = 1:NT
    % compute initial mu, Sig, S, Q, G
    phi{k} = PHI{k}(:,index);
    Hessian = alpha+phi{k}'*phi{k};
    Sig{k} = 1/Hessian;
    mu{k} = Sig{k}*PHIt(index,k);
    left = PHI{k}'*phi{k};
    S(:,k) = PHI2(:,k)-Sig{k}*left.^2;
    Q(:,k) = PHIt(:,k)-Sig{k}*PHIt(index,k)*left;
    G(:,k) = G2(:,k)-Sig{k}*PHIt(index,k)^2;
end
clear PHI2 left;
%
for count = 2:10000

    s = S; q = Q; g = G;
    Alpha = repmat(alpha,[1,NT]);
    s(index,:) = Alpha.*S(index,:)./(Alpha-S(index,:));
    q(index,:) = Alpha.*Q(index,:)./(Alpha-S(index,:));
    g(index,:) = g(index,:)+Q(index,:).^2./(Alpha-S(index,:));
    theta = NT./sum((K.*q.^2./g-s)./(s.*(s-q.^2./g)),2);

    % choice the next alpha that maximizes marginal likelihood
    ml = repmat(-inf,[M,NT]);
    ig0 = find(theta>0);
    % index for re-estimate
    [ire,foo,which] = intersect(ig0,index);
    if ~isempty(ire)
        Alpha1 = repmat(theta(ire),[1,NT]);
        Alpha0 = repmat(alpha(which),[1,NT]);
        delta = 1./Alpha1-1./Alpha0;
        X = G(ire,:).*S(ire,:)./Q(ire,:).^2;
        ml(ire,:) = (K(ire,:)-1).*log(1+S(ire,:).*delta)+K(ire,:).*log(((Alpha0+s(ire,:)).*g(ire,:)-q(ire,:).^2).*Alpha1./(((Alpha1+s(ire,:)).*g(ire,:)-q(ire,:).^2).*Alpha0));
    end
    % index for adding
    iad = setdiff(ig0,ire);
    if ~isempty(iad)
        Alpha = repmat(theta(iad),[1,NT]);
        ml(iad,:) = log(Alpha./(Alpha+s(iad,:)))-K(iad,:).*log(1-(q(iad,:).^2./g(iad,:))./(Alpha+s(iad,:)));
    end
    is0 = setdiff([1:M],ig0);
    % index for deleting
    [ide,foo,which] = intersect(is0,index);
    if ~isempty(ide)
        Alpha = repmat(alpha(which),[1,NT]);
        ml(ide,:) = -log(1-S(ide,:)./Alpha)-K(ide,:).*log(1+Q(ide,:).^2./(G(ide,:).*(Alpha-S(ide,:))));
    end
    [ML(count),idx] = max(sum(real(ml),2));

    % check if terminates?
    if count > 2 & abs(ML(count)-ML(count-1)) < (max(ML)-ML(count))*eta
        break;
    end
    % update alphas
    which = find(index==idx);
    if theta(idx) > 0
        if ~isempty(which) % re-estimate
            Alpha = theta(idx);
            delta = Alpha-alpha(which);
            for k = 1:NT
                Sigii = Sig{k}(which,which); mui = mu{k}(which); Sigi = Sig{k}(:,which);
                ki = delta/(1+Sigii*delta);
                mu{k} = mu{k}-ki*mui*Sigi;
                Sig{k} = Sig{k}-ki*Sigi*Sigi';
                comm = PHI{k}'*(phi{k}*Sigi);
                S(:,k) = S(:,k) + ki*(comm.^2);
                Q(:,k) = Q(:,k) + ki*mui*comm;
                G(:,k) = G(:,k) + ki*(Sigi'*PHIt(index,k))^2;
            end
            %
            alpha(which) = Alpha;
        else % adding
            Alpha = theta(idx);
            for k = 1:NT
                phii = PHI{k}(:,idx); Sigii = 1/(Alpha+S(idx,k)); mui = Sigii*Q(idx,k);
                comm1 = Sig{k}*(phi{k}'*phii);
                ei = phii-phi{k}*comm1;
                off = -Sigii*comm1;
                Sig{k} = [Sig{k}+Sigii*comm1*comm1', off; off', Sigii];
                mu{k} = [mu{k}-mui*comm1; mui];
                comm2 = PHI{k}'*ei;
                S(:,k) = S(:,k) - Sigii*(comm2.^2);
                Q(:,k) = Q(:,k) - mui*comm2;
                G(:,k) = G(:,k) - Sigii*(t{k}'*ei)^2;
                phi{k} = [phi{k},phii];
            end
            %
            index = [index;idx];
            alpha = [alpha;Alpha];
        end
    else
        if ~isempty(which) % deleting
            for k = 1:NT
                Sigii = Sig{k}(which,which); mui = mu{k}(which); Sigi = Sig{k}(:,which);
                Sig{k} = Sig{k}-Sigi*Sigi'/Sigii; Sig{k}(:,which) = []; Sig{k}(which,:) = [];
                mu{k}  = mu{k}-mui/Sigii*Sigi; mu{k}(which) = [];
                comm = PHI{k}'*(phi{k}*Sigi);
                S(:,k) = S(:,k) + (comm.^2)/Sigii;
                Q(:,k) = Q(:,k) + mui/Sigii*comm;
                G(:,k) = G(:,k) + (Sigi'*PHIt(index,k))^2/Sigii;
                phi{k}(:,which) = [];
            end
            %
            index(which) = [];
            alpha(which) = [];
        end
    end

end
% output
weights	= zeros(M,NT);
for k = 1:NT
    weights(index,k) = mu{k};
end
