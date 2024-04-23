function schemaball(strNames, corrMatrix, cm, xmin, xmax, cmNode, theta)
%% SCHEMABALL(strNames, corrMatrix, fontsize, positive_color_hs, negative_color_hs, theta)
%	inspired by http://mkweb.bcgsc.ca/schemaball
%	discussion at http://stackoverflow.com/questions/17038377
%
%	Draws a circular represenation of a correlation matrix.
%	Use no input arguments for a demo.
%
%INPUT ARGUMENTS
%	strNames	The names of variables of the correlation matrix
%				Format: Mx1 cell array of strings
%	corrMatrix	Correlation matrix (MxM)
%OPTIONAL:
%	fontsize	Font size of the labels along the edge
%				Default: 30/exp(M/30)
%	positive_color_hs:
%				The hue and saturation of the connection lines for
%				variables with a positive correlation.
%	negative_color_hs:
%				The hue and saturation of the connection lines for
%				variables with a negative correlation.
%   cmNode:
%               Color vector for nodes (background of text label)
%	theta		Mx1 vector of angles at which the labels and connector
%				lines must be placed. If not supplied, they are evenly
%				distributed along the whole edge of the circle.
%
% Original: https://github.com/GuntherStruyf/matlab-tools/blob/master/schemaball.m
% Adapted by Siemon Lange with tweaks by Tom Bresser

	if nargin==0 % DEMO
		strNames = {'Lorem','ipsum','dolor','sit','amet','consectetur', ...
			'adipisicing','elit','sed','do','eiusmod','tempor','incididunt',...
			'ut','labore','et','dolore','magna','aliqua'};
		
		% generate random symmetric matrix:
		a=rand(numel(strNames));
		corrMatrix=triu(a) + triu(a,1)';
		corrMatrix = 2*corrMatrix-max(corrMatrix(:));
		corrMatrix = corrMatrix.^5;
	else
		narginchk(2,7);
	end
	%% Check input arguments
	M = numel(strNames);
    if ndims(corrMatrix)~=2 || size(corrMatrix,1)~=size(corrMatrix,2) || length(corrMatrix)~=M
		error('SchemaBall:InvalidInputArguments','Invalid size of ''corrMatrix''');
    end
    
    if nargin<5
        xmin = min(corrMatrix(:));
        xmax = max(corrMatrix(:));
    elseif nargin<6
        cmNode = repmat('none', numel(strNames),1);
        theta = linspace(0,2*pi,M+1);
		theta(end)=[];
    elseif nargin<7
		theta = linspace(0,2*pi,M+1);
		theta(end)=[];
	elseif ~all(size(theta)==size(strNames))
		error('SchemaBall:InvalidInputArguments','Invalid size of ''theta''');
    end
	
	%% Configuration
	R = 1;
	Nbezier = 100;
	bezierR = 0.1;
	markerR = 0.025;
	labelR = 1.1;  
    
    fontsize = 12;
	
	%% Create figure with invisible axes, just black background
% 	figure;%('Renderer','zbuffer');
	hold on
	set(gca,'color',[1 1 1],'XTick',[],'YTick',[]);
	set(gca,'xlim',1.6*[-1 1]*R,'ylim',1.6*[-1 1]*R);       % default 2*[]*R - handles white space around circle
    %axis off
	
	%% draw diagonals
	% if you draw the brightest lines first and then the darker lines, the
	% latter will cut through the former and make it look like they have
	% holes. Therefore, sort and draw them in order (darkest first).
	idx = nonzeros(triu(reshape(1:M^2,M,M),1));
	[~,sort_idx]=sort(abs(corrMatrix(idx)));
	idx = idx(sort_idx);
	
	[Px,Py] = pol2cart(theta,R);
	P = [Px ;Py];
	
    colors = value2Color(corrMatrix(:), cm, xmin, xmax);
    colors = reshape(colors, [size(corrMatrix) 3]);
    
	for ii=idx'
		[jj,kk]=ind2sub([M M],ii);
		[P1x,P1y] = pol2cart((theta(jj)+theta(kk))/2,bezierR);
		Bt = getQuadBezier(P(:,jj),[P1x;P1y],P(:,kk), Nbezier);
% 		if corrMatrix(jj,kk)>=0
% 			clr = hsv2rgb([positive_color_hs abs(corrMatrix(jj,kk))]);
% 		else
% 			clr = hsv2rgb([negative_color_hs abs(corrMatrix(jj,kk))]);
% 		end

        %if corrMatrix(jj,kk)>0
        if corrMatrix(jj,kk)~=0   % TB - enable plotting negative numbers
            plot(Bt(:,1),Bt(:,2), ...
                'color',[squeeze(colors(jj, kk, :)); 1], ...
                'LineWidth', 7*abs(corrMatrix(jj,kk)));
        end
	end
	
	%% draw edge markers
	[Px,Py] = pol2cart(theta,R+markerR);
	% base the color of the node on the 'degree of correlation' with other
	% variables:
	corrMatrix(logical(eye(M)))=0;
	V = mean(abs(corrMatrix),2);
	V=V./max(V);
	clr = hsv2rgb([ones(M,1)*[0.585 0.5765] V(:)]);
	
	%scatter(Px,Py,20,clr);%,'filled'); % non-filled looks better imho
% 	for ii=1:M
% 		rectangle('Curvature',[1 1],'edgeColor',clr(ii,:),...
% 			'Position',[Px(ii)-markerR Py(ii)-markerR 2*markerR*[1 1]]);
% 	end
	
	%% draw labels
    resFac = 10;
    theta2 = linspace(0,2*pi,(M+1)*resFac);
    theta2(end)=[];
    
	[Px,Py] = pol2cart(theta2,labelR+0.05);
% 	for ii=1:M
% 		text(Px(ii),Py(ii),strNames{ii},'Rotation',theta(ii)*180/pi,'color',[0 0 0], ...
% 			'FontName','FixedWidth','FontSize',fontsize, 'VerticalAlignment','baseline', ...
%             'BackgroundColor', cmNode(ii,:));
%   end
    
    % draw colored circle as background  
    % split the vector into parts
    n_rows = unique(cmNode,'stable','rows');
    
    start_idx = 1;
    for i = 1:length(n_rows)
        node_index = cmNode==n_rows(i,:);
        chunk_size = sum(all(node_index,2))*resFac;
        
        end_idx = start_idx + chunk_size ; %-1
        
        if i == length(n_rows)
            Px_chunk = Px(start_idx:end);
            Py_chunk = Py(start_idx:end); 
        else
            Px_chunk = Px(start_idx:end_idx);
            Py_chunk = Py(start_idx:end_idx); 
        end


        plot(Px_chunk, Py_chunk, 'b-', 'color', n_rows(i,:), 'LineWidth', 20);
        hold on;

        start_idx = end_idx;      
    end
    
    % add labels
    [Px,Py] = pol2cart(theta,labelR);
    for ii=1:M
		text(Px(ii),Py(ii),strNames{ii},'Rotation',theta(ii)*180/pi,'color',[0 0 0], ...
			'FontName','FixedWidth','FontSize',fontsize, 'VerticalAlignment','baseline', ...
            'BackgroundColor', 'none');
    end
    
end
function Bt = getQuadBezier(p0,p1,p2,N)
	% defining Bezier Geometric Matrix B
	B = [p0(:) p1(:) p2(:)]';
	
	% Bezier Basis transformation Matrix M
	M =[1	0	0;
		-2	2	0;
		1	-2	1];
	% Calculation of Algebraic Coefficient Matrix A
	A = M*B;
	% defining t axis
	t = linspace(0,1,N)';
	T = [ones(size(t)) t t.^2];
	% calculation of value of function Bt for each value of t
	Bt = T*A;
end
