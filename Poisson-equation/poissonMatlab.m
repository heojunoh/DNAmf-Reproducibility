function [u] = poissonMatlab(a,t)
f = @(location,state)(a^2 - 2*pi^2).*exp(a*location.x).*sin(pi*location.x).*sin(pi*location.y)+2*a*pi.*exp(a.*location.x).*cos(pi*location.x).*sin(pi*location.y);

model = createpde();
geometryFromEdges(model,@squareg2);
CA = specifyCoefficients(model,'m',0,...
                               'd',0,...
                               'c',-1,...
                               'a',0,...
                               'f',f);
                           
applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',0);
generateMesh(model,'Hmax',t);

 
results = solvepde(model);

%u = mean(results.NodalSolution);
u = max(results.NodalSolution);

%[A,AE] = area(model.Mesh);
%xyEM = zeros([2, size(results.Mesh.Elements)]);
%xyEM(:) = results.Mesh.Nodes(:, results.Mesh.Elements);
%xyEc = squeeze(mean(xyEM, 2));
%sol_nodes = results.interpolateSolution(xyEc(1,:), xyEc(2,:));
%u = sum(sol_nodes'.*AE)/sum(AE);

end