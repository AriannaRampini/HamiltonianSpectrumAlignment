function plot_mask(M, N, mask)

    % plot the logical vector 'mask' on the shape 'M' and compare it with
    % the shape 'N' (it should be the red part)

    figure
    subplot(121), colormap([1 1 1])
    plot_scalar_map(N, ones(N.n,1)); hold on;
    plot_boundary_edges(N, calc_boundary_edges(N.TRIV), [0 0 1]); axis off;
    light, camlight head, lighting phong; 
    freeze_colors, view([0 90]);
    subplot(122), colormap(whitered)
    plot_scalar_map(M, mask, 0.6); axis off; %colorbar;
    light, camlight head, lighting phong, 
    view([0 90]);

end