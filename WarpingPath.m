function [m,w,c] = WarpingPath(varargin)
% (Constrained) Dynamic time warping (DTW) using a euclidean matrix and warp-path
% approach. Input vectors x & y and it returns the Euclidean matrix, m, 
% warp-path, w and cost, c.
%
% Seems useful for non-linear correlations, or coherence analysis where
% signals are almost identical but slightly out of phase.
% *Input 'demo' or 'demo2' for demo data and figures.
%
% AS2016 



% DEMO(S)
% --------------------------------
if ischar(varargin{1})
    if strcmp(varargin{1},'demo')
        
        % Generate a signal, x, & offset-by-one y
        x = randi([0 20],50,1);
        y = [0; x(1:end-1)];
        t = 1:length(x);

        % Obtain solution
        [m,w,c] = WarpingPath(x,y);
        dy      = -w*y;

        figure,
        subplot(221),plot(t,x,t, y,'--','LineWidth',3); legend({'x','y'});
        subplot(222),plot(t,x,t,dy,'--','LineWidth',3); legend({'x','-w*y'});
        subplot(2,2,[3 4]),imagesc(m);axis square
        return;
    end
    if strcmp(varargin{1},'demo2')
        
        % Generate oscillations
        x = cos(2*pi*18*(1:399)/400)'; 
        y = cos(2*pi*6*(1:399)/400)'; 
        t = 1:length(x);
        
        % Obtain solution
        [m,w,c] = WarpingPath(x,y);
        dy      = -w*y;

        figure,
        subplot(221),plot(t,x,t, y,'--','LineWidth',3); legend({'x','y'});
        subplot(222),plot(t,x,t,dy,'--','LineWidth',3); legend({'x','-w*y'});
        subplot(2,2,[3 4]),imagesc(m);axis square
        return;
    end
end



% The Euclidean function
%--------------------------------
try x = varargin{1}; end
try y = varargin{2}; catch y = x; end

% ensure vectors
x = x(:); xl = length(x);
y = y(:); yl = length(y);

% euclidean [warping] matrix
for i = 1:xl
    for j  = 1:yl
        x0 = x(i);
        y0 = y(j);     
        E(i,j) = sum ( (x0-y0).^2 );  
    end
end

% check whether to find path or give up
if nargout == 1; m = E; return; end



% The Path function
%--------------------------------

% expansion point for search
P   = E(1,1); 
C   = 0;
V   = 0;
Coo = [1,1];

xl1   = xl;
for i = 2:xl*.5*xl
    j = i;
    
    % update current co-ordinate
     fx = Coo(1); fy = Coo(2);
    
    % check for matrix completion
    if size(P,2) == xl1 || size(P,1) == xl1
        if     size(P,1) > size(P,2); P(xl,xl) = E(end,end); break
        elseif size(P,1) < size(P,2); P(xl,xl) = E(end,end); break
        end
    end

    % what are the 3 adjacent options:
    choice{1} = E(fx+1,fy  ); % dx/ y
    choice{2} = E(fx+1,fy+1); % dx/dy
    choice{3} = E(fx  ,fy+1); %  x/dy
    
    % ensure it can't stray too far from diag [by reducing options]
    if abs(fx-fy) >= round(1/16*xl1)
        if     fx > fy; choice{1} = Inf;
        elseif fy > fx; choice{3} = Inf;
        end
    end
    
    % select minimum distance & log
    [v,op] = min(cat(1,choice{:}));
    
    % if same exact value
    if v == 0 
       v = -1;
    end
    
    C(i)   = op; 
    V(i)   = v; 
    
    % cost function
    cost(i) = abs(Coo(1)-Coo(2))/i;

    % new coordinates
    if     op == 1; Coo(1) = Coo(1) + 1;
    elseif op == 2; Coo = Coo + 1;
    elseif op == 3; Coo(2) = Coo(2) + 1;
    end
    
    % record value in matrix    
    if     op == 1;           P(fx+1,fy  ) = v; 
    elseif op == 2 && i ~= 2
        if isempty(P(fx,fy)); P(fx,fy)     = v; 
        else                  P(fx+1,fy+1) = v;
        end
    elseif op == 2 && i == 2; P(fx+1,fy+1) = v;
    elseif op == 3;           P(fx  ,fy+1) = v; 
    end

    
end


m = E;
w = P;
c = cost;
