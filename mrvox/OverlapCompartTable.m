function OverlapMatrice = OverlapCompartTable(Centers,Radius)


C = Centers;
R = Radius;
n = size(C,1);

dist = ...
    (( C(:,1)*ones(1,n) - ones(n,1)*C(:,1)' ) .^2 + ...
     ( C(:,2)*ones(1,n) - ones(n,1)*C(:,2)' ) .^2 ).^(1/2);
 
dist(dist==0) = Inf;

OverlapMatrice = dist < ones(n,1)*R(1:n) + (ones(n,1)*R(1:n))';
