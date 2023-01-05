% Relative error for the unknown u
figure
semilogy(1:R,ef_u_pod,'-b',1:R,ef_u_deim,'-dm',1:R,ef_u_podc,'-.r',1:R,ef_u_deimc,'-.dg','LineWidth',2,'MarkerSize',2)
xlabel('r')
legend('POD','POD-DEIM','PODc','POD-DEIMc')
title('$\cal{E}({\bf{u}}$,r)','Interpreter','latex')
set(gca,'FontSize',16)

% Relative error for the unknown v
figure
semilogy(1:R,ef_v_pod,'-b',1:R,ef_v_deim,'-dm',1:R,ef_v_podc,'-.r',1:R,ef_v_deimc,'-.dg','LineWidth',2,'MarkerSize',2)
xlabel('r')
legend('POD','POD-DEIM','PODc','POD-DEIMc')
title('$\cal{E}({\bf{v}}$,r)','Interpreter','latex')
set(gca,'FontSize',16)

% Compare CPU time (s)
figure
semilogy(1:R,time_pod,'-b',1:R,time_deim,'-dm',1:R,time_podc,'-.r',1:R,time_deimc,'-.dg','LineWidth',2,'MarkerSize',2)
xlabel('r')
ylabel('CPU time (s)')
legend('POD','POD-DEIM','PODc','POD-DEIMc')
set(gca,'FontSize',16)    