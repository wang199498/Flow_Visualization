
%posuchen Wang
%v3   3D
%generate paticles in box #=2000
%trace their motion bpos XXX
%delete some of them and random add some
%[73, 178, 192]./255

clear
close all
clc
%%%Programmed Flow List
library={'Lorenzo'};
% flow_name=library{randi(length(library))};
flow_name=library{1};

%%% visualization parameters
moviebox_size=40;
n_particles=250;
max_time_steps=500;
dt=0.005;
delete_ratio=0.01;
flow_memory=200;
delete_num=round(n_particles*delete_ratio);


%%%Particle position initializition
pos_stock=cell(max_time_steps,1);
pos0=rand(n_particles,3)*moviebox_size*2-moviebox_size;
pos0(:,3)=rand(n_particles,1)*50;
pos_stock{1}=pos0;
%%%Compute all positions along time
for time_step=2:max_time_steps
    %%% define the vector field
    x=pos0(:,1);
    y=pos0(:,2);
    z=pos0(:,3);
    if strcmp(flow_name,'Lorenzo')
        vx=10*(y-x);
        vy=x.*(28-z)-y;
        vz=x.*y-8/3*z;
    elseif strcmp(flow_name,'Sink')
        vx=-x./((x).^2+y.^2);
        vy=-y./(x.^2+y.^2);
    else
        vx=0.1*sin(y);
        vy=px;
    end

    %%% forward Euler propagation
    pos(:,1)=pos0(:,1)+dt*vx;
    pos(:,2)=pos0(:,2)+dt*vy;
    pos(:,3)=pos0(:,3)+dt*vz;
    %%% patrial refresh data
    pos_cut=pos(delete_num+1:end,:);
    pos_add=rand(delete_num,3)*moviebox_size*2-moviebox_size;
    pos=[pos_cut; pos_add];
    %%% write pos data to stock
    pos_stock{time_step}=pos;
    %%% update for iteration
    pos0=pos;
end

%%% color play
c=cell(max_time_steps,1);
for time_step=1:max_time_steps
    c{time_step}=zeros(n_particles,3);
    for i=1:n_particles
        % %%uniform ditribution color
%         c{time_step}(i,:)=[73, 178, 192]./255;
        %%cicluar ditribution color (apeture)
                r=sqrt((pos_stock{time_step}(i,1))^2+(pos_stock{time_step}(i,2))^2+((pos_stock{time_step}(i,3)-25)^2));
        

                    c{time_step}(i,:)=r/20*[73, 178, 192]./255+(1-r/20)*[1 1 1];
               
        % % %%binary ditribution color
        %         if (pos_stock{time_step}(i,1))<0
        %             c{time_step}(i,:)=[0, 0, 0];
        %         else
        %             c{time_step}(i,:)=[1, 1, 1];
        %         end
        % %%binary gradient ditribution color
%         dis=(pos_stock{time_step}(i,1));
%         dis=(dis+1)/2;
%         c{time_step}(i,:)=dis*[73, 178, 192]./255+(1-dis)*[1 1 1];


    end
end
%%
clear pp
% Plot as movie
f=figure;
f.Position = [400 0 1000 1000];
xlim([-moviebox_size,moviebox_size])
ylim([-moviebox_size,moviebox_size])
zlim([0 50])
% xticklabels({''})
% yticklabels({''})
% zticklabels({''})
% axis square
title(strcat({'Flow Type='},{' '},flow_name))
set(gca,'Color','#242C39')    %sky in the dark
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',30,'FontWeight','bold')
% set(f.axes1,'XLim',[-10 10],'posLim',[-10 10]);
hold all
for time_step=1:length(pos_stock)
     pp(time_step)=scatter3(pos_stock{time_step}(:,1),pos_stock{time_step}(:,2),pos_stock{time_step}(:,3),0.2,c{time_step});
    view(32,55)
%     view(32*time_step/max_time_steps,55)
    if length(pp)>flow_memory
        delete(pp(time_step-flow_memory))
    end
    
    pause(0.01)

end
