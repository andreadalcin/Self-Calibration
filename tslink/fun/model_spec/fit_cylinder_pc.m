function [m] = fit_cylinder_pc(X,epsi)
%FIT_CYLINDER_PC use the pc library to fit cylinder
 pc = pointCloud(X(1:3,:)','Normal',X(4:6,:)');
 model = pcfitcylinder(pc,epsi);
 finiteCyinderlPar = model.Parameters(:);
 m = convertFiniteCylinderToParam(finiteCyinderlPar);
 
end

