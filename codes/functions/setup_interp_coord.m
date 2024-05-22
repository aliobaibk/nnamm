function [X, Y, X_q_norm, Y_q_norm, X_q, Y_q, m_X, m_Y, s_X, s_Y] = ...
    setup_interp_coord(X, Y, n_X_q, n_Y_q)
X = X(:);
m_X = mean(X);
s_X = std(X);
X = (X - m_X) / s_X;

Y = Y(:);
m_Y = mean(Y);
s_Y = std(Y);
Y = (Y - m_Y) / s_Y;

[X_q_norm, Y_q_norm] = meshgrid(linspace(min(X), max(X), n_X_q), linspace( ...
    min(Y), max(Y), n_Y_q));
X_q = X_q_norm * s_X + m_X;
Y_q = Y_q_norm * s_Y + m_Y;
end
