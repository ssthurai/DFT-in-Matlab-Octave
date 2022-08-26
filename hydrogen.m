g = 100; g3 = g^3; %grid points and grid dimensions

p = linspace(-8,8,g); %defining grid for one-dimension

%Now defining grid for 3-dimensions.
[X,Y,Z] = meshgrid(p,p,p);

%Now grid spacing  
h = p(2) - p(1);

%making it a one-dimensional
X=X(:); Y=Y(:); Z=Z(:);

%distance of electron

R = sqrt(X.^2+Y.^2+Z.^2);

%Now the external potential
Vext = -1./R;

%Now we need to calculate the kinetic energy. 

e=ones(g,1);  %helper identity matrix

L = spdiags([e -2*e e], -1:1, g, g)/h^2;    #Laplacian operator

I = speye(g);  %helper sparse identity matrix

L3 = kron(kron(L,I),I) + kron(kron(I,L),I) + kron(kron(I,I),L);

E = eigs(-0.5*L3+spdiags(Vext, 0, g3, g3), 1, 'sa');

disp(['Total energy for H atom ' num2str(E*27.21,6) ' eV.']);


 
