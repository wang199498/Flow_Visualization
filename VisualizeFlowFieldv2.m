
%posuchen Wang
%v2
%generate paticles in box #=2000
%trace their motion bpos XXX
%delete some of them and random add some
%[73, 178, 192]./255

clear
close all
clc
%%%Programmed Flow List
library={'Sink',...
    'Source',...
    'Rotating Source',...
    'Doublet',...
    'Quadrupole',...
    'Star Trace'...
    'Multi Source',...
    'Lab1',...
    'Lab2',...
    'Lab3',...
    'Lab4',...
    '?'};
% flow_name=library{randi(length(library))};
flow_name=library{end-1};

%%% visualization parameters
% colorlist={'#7EAA8F','#A6C7B3','#35414A'}
% colorlist={'r','g','b'};
moviebox_size=1;
n_particles=2000;
max_time_steps=1000;
% dt=0.005;
dt=0.09;
delete_ratio=0.05;
flow_memory=100;
delete_num=round(n_particles*delete_ratio);


%%%Particle position initializition
pos_stock=cell(max_time_steps,1);
pos0=rand(n_particles,2)*moviebox_size*2-moviebox_size;
pos_stock{1}=pos0;
%%%Compute all positions along time
for time_step=2:max_time_steps
    %%% define the vector field
    x=pos0(:,1);
    y=pos0(:,2);

    if strcmp(flow_name,'Source')
        vx=x./((x).^2+y.^2);
        vy=y./(x.^2+y.^2);
    elseif strcmp(flow_name,'Sink')
        vx=-x./((x).^2+y.^2);
        vy=-y./(x.^2+y.^2);
    elseif strcmp(flow_name,'Rotating Source')
        vx=-y./((x).^2+y.^2)+x./((x).^2+y.^2)-sin(time_step*0.05);
        vy=x./(x.^2+y.^2)+y./(x.^2+y.^2)-cos(time_step*0.05);
    elseif strcmp(flow_name,'Doublet')
        vx=(y.^2-x.^2)./(x.^2+y.^2);
        vy=-(2.*x.*y)./(x.^2+y.^2);
    elseif strcmp(flow_name, 'Quadrupole')
        r=sqrt(x.^2+y.^2);
        ur=-2*(1-2*x.^2./(r.^2))./(r.^3)*0.1;
        ut=-(2.*x.*y./(r.^2))./(r.^3)*0.1;
        vx=ur.*y./r-ut.*x./r;
        vy=ur.*x./r+ut.*y./r;
    elseif strcmp(flow_name,'Star Trace')
        vx=-y;
        vy=x;
    elseif strcmp(flow_name,'Multi Source')
        [sx1,sy1]=source_part(x+0.5,y+0.5);
        [sx2,sy2]=source_part(x+0.5,y-0.5);
        [sx3,sy3]=source_part(x-0.5,y+0.5);
        [sx4,sy4]=source_part(x-0.5,y-0.5);
        vx=sx1+sx2+sx3+sx4;
        vy=sy1+sy2+sy3+sy4;
    elseif strcmp(flow_name,'Lab1')
        vx=sin(5*y );
        vy=cos(5*x );
    elseif strcmp(flow_name,'Lab2')
        d = sqrt(x.^2+y.^2);
        vx = 2*cos(0.4*pi*time_step*0.01./d);
        vy = 2*sin(0.4*pi*time_step*0.01./d);
    elseif strcmp(flow_name,'Lab3')
        range1=norm(y);
        range2=norm(x+6*sign(x));
        vx = -2.6*y;
        vy = 2.6*(x+6*sign(x));
            elseif strcmp(flow_name,'Lab4')

        vx =  y;
%         vy = 0.5*(1-x.^2).*y-x;

vy= -y+x-x.^3;
    else
        vx=0.1*sin(y);
        vy=px;
    end

    %%% forward Euler propagation
    pos(:,1)=pos0(:,1)+dt*vx;
    pos(:,2)=pos0(:,2)+dt*vy;
    %%% patrial refresh data
    pos_cut=pos(delete_num+1:end,:);
    pos_add=rand(delete_num,2)*moviebox_size*2-moviebox_size;
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
%                 r=sqrt((pos_stock{time_step}(i,1))^2+(pos_stock{time_step}(i,2))^2);
%         
%                 if r>0.5
%                     c{time_step}(i,:)=[0, 0, 0];
%                 else
%                     c{time_step}(i,:)=r/0.5*[73, 178, 192]./255+(1-r/0.5)*[1 1 1];
%                 end
        % % %%binary ditribution color
        %         if (pos_stock{time_step}(i,1))<0
        %             c{time_step}(i,:)=[0, 0, 0];
        %         else
        %             c{time_step}(i,:)=[1, 1, 1];
        %         end
        % %%binary gradient ditribution color
        dis=(pos_stock{time_step}(i,1));
        dis=(dis+1)/2;
        c{time_step}(i,:)=dis*[73, 178, 192]./255+(1-dis)*[1 1 1];


    end
end
%%
clear pp
% Plot as movie
f=figure;
f.Position = [400 0 1000 1000];
xlim([-moviebox_size,moviebox_size])
ylim([-moviebox_size,moviebox_size])
xticklabels({''})
yticklabels({''})
% axis square
title(strcat({'Flow Type='},{' '},flow_name))
set(gca,'Color','#242C39')    %sky in the dark
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',30,'FontWeight','bold')
% set(f.axes1,'XLim',[-10 10],'posLim',[-10 10]);
hold all
for time_step=1:length(pos_stock)
    p=round(rand*2+1);
    %         pp(time_step)=plot(pos_stock{time_step}(:,1),pos_stock{time_step}(:,2),'o','markersize',0.2,'color',[73, 178, 192]./255);
    %     pp(time_step)=plot(pos_stock{time_step}(:,1),pos_stock{time_step}(:,2),'o','markersize',0.2,'color',rand(3,1));
    % c=1:1:length(pos_stock{time_step}(:,1))   ;



    pp(time_step)=scatter(pos_stock{time_step}(:,1),pos_stock{time_step}(:,2),0.2,c{time_step});

    % colormap(gca,'winter')


    if length(pp)>flow_memory
        delete(pp(time_step-flow_memory))
    end
    pause(0.01)

end



function [sx,sy]=source_part(x,y)
sx=x./((x).^2+y.^2);
sy=y./(x.^2+y.^2);
end

function f=smoothstep(edge1,edge2,x)
if x<=edge1
    f=0;
elseif x>=edge2
    f=1;
else
    f=3*x.^2-2*x^3;
end
end
