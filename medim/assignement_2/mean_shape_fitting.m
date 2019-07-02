%val_models = models(:,1:5);
%ms_x = mean_shape(1:14,:);
%ms_y = mean_shape(1:14,:);
results = nan(size(models, 2), 2);
c = 1;
ms_x = mean_shape(1:14,:);
ms_y = mean_shape(15:end,:);
for vm = models
    vm_x = real(vm);
    vm_y = imag(vm);
    
    [~,rot_mean] = procrustes([vm_x, vm_y],[ms_x, ms_y]);
    rot_mean_ps = polyshape(rot_mean(:,1), rot_mean(:,2));
    val_ps = polyshape(vm_x, vm_y);
    
    results(c, 1) = area(intersect(val_ps, rot_mean_ps)) / area(val_ps);
    results(c, 2) = area(xor(val_ps, rot_mean_ps)) / area(val_ps);
    c = c+1;
end