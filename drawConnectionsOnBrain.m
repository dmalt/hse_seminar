function drawConnectionsOnBrain(ConInds, GridXYZ, Ctx)
% -------------------------------------------------------
% Given indices of connected grid locations, grid nodes
% coordinates and cortex surface, draw connections on
% brain
% -------------------------------------------------------
% FORMAT:
%   drawConnectionsOnBrain(ConInds, GridXYZ, iCol, Ctx) 
% INPUTS:
%   ConInds        - {nConnections x 2} matrix of indices
%                    of connected grid nodes
%   GridXYZ		   - {nSources x 3} matrix of coordinates
%                    of grid nodes
%   iCol           - int; color number
%   Ctx            - structure; triangular surface of cortex;
%                    usually generated in brainstorm
% OUTPUTS:
% _______________________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru
	grey =	[242, 237, 217] / 256;
    color = [241, 169, 78] / 256;
	figure;
	trisurf(Ctx.Faces,Ctx.Vertices(:,1),Ctx.Vertices(:,2),Ctx.Vertices(:,3),...
					'FaceColor', grey, 'EdgeColor','none','FaceAlpha', 0.2); 
	lighting phong;
	camlight right;

	set(gcf,'color','k');
	axis equal;
	view(-90,90);
	axis off;
	grid off;

	hold on;
	drawset(ConInds, GridXYZ, color);
end
