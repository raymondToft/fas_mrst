function h = plotCircle(r, varargin)
    center = [0, 0];
    if numel(varargin) && isnumeric(varargin{1})
        center = varargin{1};
        varargin = varargin(2:end);
    end
    if numel(r) == 1
        r = [r, r];
    end
    a = 0:pi/50:2*pi;
    h = plot(center(1) + r(:, 1).*cos(a), center(2) + r(:, 2).*sin(a), varargin{:});
end