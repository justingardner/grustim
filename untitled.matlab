pt = {'Brillantes','Nguyen','Sih','Valandani'};
figure;
for p = 1:length(pt)
	curdat = eval(pt{p});
	ind = find(strcmp(curdat(:,1),'INC'));
	resp = curdat(ind,3);
	nTrials = length(resp);
	nb = sum(strcmp(resp,'b'));
	nd = sum(strcmp(resp,'d'));
	ng = sum(strcmp(resp,'g'));

	pb(p) = nb/nTrials;
	pd(p) = nd/nTrials;
	pg(p) = ng/nTrials;

end

scatter3(pb,pg,pd,150,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[.1 .5 .75]);
axis([0 1 0 1 0 1]);
xlabel('Auditory Dominance (% "ba")','FontSize',14);
zlabel('Bimodal Fusion (% "da")', 'FontSize',14);
ylabel('Visual Dominance (% "ga")', 'FontSize',14);

