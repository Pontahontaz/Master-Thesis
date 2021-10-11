function g = odampedlineshape(t_ls,lam,tau,T)
%LINESHAPE Generates the lineshape function
% g = lineshape(t_ls,lam,tau,HR_vibr,tau_vibr,nu_vibr,T)
%   Butkus J. Chem. Phys. 137:044513 (2012)

% conversion factors
cm2w = 2 * pi / 33356.4095; % [PHz rad cm]

% lineshape overdamped parameters
gam = 2 * pi ./ (tau*1000); % [PHz rad]
lam = lam * cm2w; % [PHz rad]


% additional constants
kB = 1.3806488e-23; % [J K^{-1}]
h = 6.62606957e-34 * 1e15; % [J fs]
beta = 1 / (kB * T); % [J^{-1}]
hbar = h / (2 * pi); % [J rad^{-1} fs]
beta = beta * hbar; % [fs rad^{-1}]

% overdamped part
g_over = zeros( size(t_ls) );
for i = 1:numel(gam)
        g_over = g_over + lam(i) / gam(i) * (2 / (beta * gam(i)) - 1i) *...
            (exp(- gam(i) * t_ls) + gam(i) * t_ls - 1);
%         g_over = g_over + 0.5*(gam(i)*t_ls).^2 + gam(i)*t_ls; % effective lineshape
end

g.x = permute(g_over,[1 3 2]);
g.t = t_ls;

end