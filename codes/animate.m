function animate(phi, omega, varargin)

% ------------ Parse ------------
assert(size(phi,1)==10, 'phi must be 10xM with 10 rows.');
M = size(phi,2);  assert(M>=1 && M<=3, 'Provide 1–3 modes (columns).');

p = inputParser;
p.addParameter('Scale', [], @(x) isempty(x) || isscalar(x) || (isvector(x) && numel(x)==M));
p.addParameter('Periods', 5, @(x) isscalar(x) && x>0);
p.addParameter('FPS', 60, @(x) isscalar(x) && x>0);
p.addParameter('BlackAlpha', 0.72, @(x) isscalar(x) && x>=0 && x<=1);
p.addParameter('LineAlpha', [], @(x) isempty(x) || (isscalar(x) && x>=0 && x<=1));
p.addParameter('BlockAlpha', [], @(x) isempty(x) || (isscalar(x) && x>=0 && x<=1));
p.addParameter('FigurePosition', [100 100 1500 520], @(v) isnumeric(v) && numel(v)==4);
p.addParameter('PaperSize', [22 12], @(v) isnumeric(v) && numel(v)==2);
p.addParameter('Renderer', 'opengl', @(s) ischar(s) || isstring(s));
p.addParameter('ModeTitles', {}, @(c) iscellstr(c) || isstring(c));
p.parse(varargin{:});
opt = p.Results;

% normalize omega to row vector 1xM
if isscalar(omega), omega = repmat(omega,1,M);
else, omega = omega(:).';  assert(numel(omega)==M,'omega must be scalar or length M.');
end

% Resolve alphas
aLine  = opt.BlackAlpha;  if ~isempty(opt.LineAlpha),  aLine  = opt.LineAlpha;  end
aBlock = opt.BlackAlpha;  if ~isempty(opt.BlockAlpha), aBlock = opt.BlockAlpha; end

% Titles
if isempty(opt.ModeTitles)
    opt.ModeTitles = repmat({'Normal mode representation: eigenvector'},1,M);
else
    opt.ModeTitles = cellstr(opt.ModeTitles);
    if numel(opt.ModeTitles) < M
        opt.ModeTitles(end+1:M) = opt.ModeTitles(end);
    end
end

% ------------ Geometry (shared) ------------
xR =  1;  xL = -1;
y1=0; y2=2; y3=4; y4=5; y5=6; y6=7;
y7=3; y8=4; yNode=5; y9=6; y10=7;
x0 = [xR xR xR xR xR xR  xL xL xL xL];
y0 = [y1 y2 y3 y4 y5 y6  y7 y8 y9 y10];
yNode0 = yNode;
cLine  = [0 0 0];
cGhost = [0.7 0 0];

spring_coords = @(x,yA,yB) deal( ...
    (x + 0.18 * (-1).^(1:15)).*[1  ones(1,13) 1], ... % will fix ends after
    linspace(yA,yB,15) );
% (fix ends in a helper)
    function [xx,yy]=spr(x,yA,yB)
        [xx,yy] = spring_coords(x,yA,yB);
        xx(1)=x; xx(end)=x;
    end

% ------------ Figure scaffold ------------
fig = gcf;
set(fig,'Units','pixels','Position',opt.FigurePosition,'Color','w', ...
         'Renderer',opt.Renderer,'InvertHardcopy','off');
set(fig,'PaperUnits','centimeters','PaperPosition',[0 0 opt.PaperSize], ...
         'PaperSize',opt.PaperSize);
clf(fig);
tl = tiledlayout(fig,1,M,'TileSpacing','compact','Padding','compact');

% Per-mode scale (auto if empty or scalar expand)
if isempty(opt.Scale)
    scale = zeros(1,M);
    for k=1:M
        scale(k) = 0.45 / max(max(abs(phi(:,k))), eps);
    end
elseif isscalar(opt.Scale)
    scale = repmat(opt.Scale,1,M);
else
    scale = opt.Scale(:).';
end

% ------------ Build each panel once, keep handles ------------
Panels = struct([]);
for k = 1:M
    ax = nexttile(tl,k); hold(ax,'on'); axis(ax,'equal'); box(ax,'on'); grid(ax,'on');
    title(ax, opt.ModeTitles{k});
    xlabel(ax,'Schematic position'); ylabel(ax,'Final vertical position');
    xlim(ax,[-2 2]); ylim(ax,[min(y0)-0.9, max(y0)+0.9]); axis(ax,'manual');

    % Ghost springs
    [xx,yy]=spr(xR,y1,y2);  plot(ax,xx,yy,'-','Color',[cGhost 0.25]);
    [xx,yy]=spr(xR,y2,y3);  plot(ax,xx,yy,'-','Color',[cGhost 0.25]);
    [xx,yy]=spr(xR,y3,y4);  plot(ax,xx,yy,'-','Color',[cGhost 0.25]);
    [xx,yy]=spr(xR,y4,y5);  plot(ax,xx,yy,'-','Color',[cGhost 0.25]);
    [xx,yy]=spr(xR,y5,y6);  plot(ax,xx,yy,'-','Color',[cGhost 0.25]);
    [xx,yy]=spr(xL,y7,y8);          plot(ax,xx,yy,'-','Color',[cGhost 0.25]);
    [xx,yy]=spr(xL,y8,yNode0);      plot(ax,xx,yy,'-','Color',[cGhost 0.25]);
    [xx,yy]=spr(xL,yNode0,y9);      plot(ax,xx,yy,'-','Color',[cGhost 0.25]);
    [xx,yy]=spr(xL,y9,y10);         plot(ax,xx,yy,'-','Color',[cGhost 0.25]);
    plot(ax,[xR xL],[y4 yNode0],'-','Color',[cGhost 0.25],'LineWidth',1.0); % ghost link

    % Deformed masses + labels (placeholders at undeformed; updated in loop)
    hMass = gobjects(10,1); hLab = gobjects(10,1);
    for i=1:10
        hMass(i)=rectangle(ax,'Position',[x0(i)-0.35, y0(i)-0.2, 0.7, 0.4], ...
            'Curvature',0.05,'EdgeColor',cLine,'FaceColor',[cLine aBlock],'LineWidth',1.2);
        hLab(i)=text(ax,x0(i)+0.45, y0(i), sprintf('x_{%d}',i),'FontSize',11, ...
            'Color',cLine,'HorizontalAlignment','left','VerticalAlignment','middle');
    end

    % Deformed springs (handles to update)
    hSpr = gobjects(9,1);
    [xx,yy]=spr(xR,y1,y2);  hSpr(1)=plot(ax,xx,yy,'-','Color',[cLine aLine]);
    [xx,yy]=spr(xR,y2,y3);  hSpr(2)=plot(ax,xx,yy,'-','Color',[cLine aLine]);
    [xx,yy]=spr(xR,y3,y4);  hSpr(3)=plot(ax,xx,yy,'-','Color',[cLine aLine]);
    [xx,yy]=spr(xR,y4,y5);  hSpr(4)=plot(ax,xx,yy,'-','Color',[cLine aLine]);
    [xx,yy]=spr(xR,y5,y6);  hSpr(5)=plot(ax,xx,yy,'-','Color',[cLine aLine]);
    [xx,yy]=spr(xL,y7,y8);          hSpr(6)=plot(ax,xx,yy,'-','Color',[cLine aLine]);
    [xx,yy]=spr(xL,y8,yNode0);      hSpr(7)=plot(ax,xx,yy,'-','Color',[cLine aLine]);
    [xx,yy]=spr(xL,yNode0,y9);      hSpr(8)=plot(ax,xx,yy,'-','Color',[cLine aLine]);
    [xx,yy]=spr(xL,y9,y10);         hSpr(9)=plot(ax,xx,yy,'-','Color',[cLine aLine]);

    % Rigid link (handle)
    hLink = plot(ax,[xR xL],[y4 yNode0],'-','Color',[cLine aLine],'LineWidth',1.2);

    % store
    Panels(k).ax = ax;
    Panels(k).hMass = hMass;
    Panels(k).hLab  = hLab;
    Panels(k).hSpr  = hSpr;
    Panels(k).hLink = hLink;
end

% ------------ Time base ------------
% Use the largest period among modes to compute a total time that covers
% 'Periods' full cycles of the **slowest** mode (so all run long enough).
omegaSafe = omega; omegaSafe(omegaSafe<=eps) = 2*pi; % avoid zero
Tslow = 2*pi / min(omegaSafe);
tEnd  = opt.Periods * Tslow;
dt    = 1/opt.FPS;
t     = 0:dt:tEnd;

% ------------ Animate all panels together ------------
for kFrame = 1:numel(t)
    tk = t(kFrame);
    for m = 1:M
        s = sin(omega(m) * tk);  % phase for panel m
        y     = y0 + scale(m) * (phi(:,m).'* s);
        yNode = yNode0 + scale(m) * phi(4,m) * s;

        % masses + labels
        for i = 1:10
            set(Panels(m).hMass(i),'Position',[x0(i)-0.35, y(i)-0.2, 0.7, 0.4]);
            set(Panels(m).hLab(i), 'Position',[x0(i)+0.45, y(i), 0]);
        end
        % springs
        [xx,yy]=spr(xR,y(1),y(2));  set(Panels(m).hSpr(1),'XData',xx,'YData',yy);
        [xx,yy]=spr(xR,y(2),y(3));  set(Panels(m).hSpr(2),'XData',xx,'YData',yy);
        [xx,yy]=spr(xR,y(3),y(4));  set(Panels(m).hSpr(3),'XData',xx,'YData',yy);
        [xx,yy]=spr(xR,y(4),y(5));  set(Panels(m).hSpr(4),'XData',xx,'YData',yy);
        [xx,yy]=spr(xR,y(5),y(6));  set(Panels(m).hSpr(5),'XData',xx,'YData',yy);
        [xx,yy]=spr(xL,y(7),y(8));  set(Panels(m).hSpr(6),'XData',xx,'YData',yy);
        [xx,yy]=spr(xL,y(8),yNode); set(Panels(m).hSpr(7),'XData',xx,'YData',yy);
        [xx,yy]=spr(xL,yNode,y(9)); set(Panels(m).hSpr(8),'XData',xx,'YData',yy);
        [xx,yy]=spr(xL,y(9),y(10)); set(Panels(m).hSpr(9),'XData',xx,'YData',yy);
        % link
        set(Panels(m).hLink,'XData',[xR xL],'YData',[y(4) yNode]);
    end
    drawnow;
end
end
