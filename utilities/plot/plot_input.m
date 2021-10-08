function plot_input(M, N)

    % plot the template 'M' and the partial shape 'N' 

    figure
    subplot(121), colormap([1 1 1])
    plot_scalar_map(N, ones(N.n,1)); hold on;
    plot_boundary_edges(N, calc_boundary_edges(N.TRIV), [0 0 1]); axis off;
    light, camlight head, lighting phong; 
    freeze_colors;  
    view([-37 19]);
    title('Partial shape N')
    
    subplot(122), colormap([1 1 1])
    plot_scalar_map(M, ones(M.n,1)); 
    axis off; 
    light, camlight head, lighting phong, 
    view([-37 19]);
    title('Template M')

end