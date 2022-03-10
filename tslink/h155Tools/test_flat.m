B = rand(3,2);
origin = rand(3,1);
%%
n = 100;
for i = 1:n
    X(:,i) = origin + rand(1).*B(:,1) + rand(1).*B(:,2)
end

flat = fit_flat(X,2);

B = flat.basis;
X_recon = ((B*B')*(X-flat.origin))+ flat.origin;



figure; 
scatter3(X(1,:),X(2,:),X(3,:),'r*');
hold on
scatter3(X_recon(1,:),X_recon(2,:),X_recon(3,:),'bo');


d = res_flat(X,flat); 


%% test the distance from a point
B = eye(3,2);
origin = zeros(3,1);

n = 100;
for i = 1:n
    X(:,i) = origin + rand(1).*B(:,1) + rand(1).*B(:,2)
end

flat = fit_flat(X,2);

B = flat.basis;
X_recon = ((B*B')*(X-flat.origin))+ flat.origin;



figure; 
scatter3(X(1,:),X(2,:),X(3,:),'r*');
hold on
scatter3(X_recon(1,:),X_recon(2,:),X_recon(3,:),'bo');

P = [0;0; 100];


d1 = res_flat(P,flat);
d2 = res_flat(P,flatStruct2Vect(flat));
