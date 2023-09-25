
set(gca,'LineWidth',1.5,'FontSize',25,'FontWeight','normal','FontName','Times')
set(gcf,'Position',[1 1 round(1000) round(1000)])
hline=findobj(gcf,'type','line');
set(hline,'Linewidth',2)
set(get(gca,'xlabel'),'String','No of high fidelity points','FontSize',28','FontWeight','bold','FontName','Times','Interpreter','tex')
%legend({'Actual system','MF-HPCFE','LF-HPCFE','HF-HPCFE'},'FontSize',20,'Location','southeast')
%legend('FontSize',20,'Location','southeast')
%legend boxoff
title("")
set(get(gca,'ylabel'),'String','Root mean square error','FontSize',28','FontWeight','bold','FontName','Times','Interpreter','tex')
set(gcf,'color','w')
box on

%export_fig(gcf,'50_16_twofidelity.pdf','-pdf','-r300')
% set(gca,'LineWidth',2,'FontSize',24,'FontWeight','normal','FontName','Times')
% set(gcf,'Position',[1 1 round(sz1) round(sz2)])
% st(get(gca,'xlabel'),'String','z_1','FontSize',32','FontWeight','bold','FontName','Times','Interpreter','tex')
% set(get(gca,'ylabel'),'String','z_2','FontSize',32','FontWeight','bold','FontName','Times','Interpreter','tex')