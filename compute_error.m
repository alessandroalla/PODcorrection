function [error] = compute_error(u_app,u_ex)

error = norm(u_app-u_ex,'fro')/norm(u_ex,'fro');

end

