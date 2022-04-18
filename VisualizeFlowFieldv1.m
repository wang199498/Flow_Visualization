
%Yuchen Wang 
% v1
%generate paticles in box #=2000
%trace their motion by ode45
%draw back, particles will flow out of box

clear
close all
clc
moviebox_size=10;
n_particles=2000;
color_list=rand(n_particles,3);
y=cell(n_particles,1);
step_stock=[];
opts = odeset('RelTol',1e-3,'AbsTol',1e-8);
for i=1:n_particles
    y0=[rand*moviebox_size*2-moviebox_size,rand*moviebox_size*2-moviebox_size];
    [t,y{i}] = ode45(@field,[0 30],y0,opts);
    step_stock=[step_stock; length(y{i})];
end
movie_steps=min(step_stock);
x_rectify=zeros(movie_steps, n_particles);
y_rectify=zeros(movie_steps, n_particles);
for i=1:n_particles
    
    x_rectify(:,i)=y{i}(1:movie_steps,1);
    y_rectify(:,i)=y{i}(1:movie_steps,2);
end
f = figure;
f.Position = [400 50 700 700];
title('Flow Field Visualization')
xlim([-moviebox_size moviebox_size])
ylim([-moviebox_size moviebox_size])
xticklabels({});
yticklabels({});
axis square
set(gca,'Color','k')
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',30,'FontWeight','bold')
hold all
for i=1:movie_steps

    pp(i)=scatter(x_rectify(i,:),y_rectify(i,:),1);
%     pp(i).MarkerEdgeAlpha=(movie_steps-i)/movie_steps;
    pp(i).MarkerEdgeColor='w';
    progress=floor(i/movie_steps*100);
    dyn_label=strcat('Progress=',num2str(progress),'%');
    xlabel(dyn_label)
    treshold_length=5;
    if length(pp)>treshold_length
        delete(pp(i-treshold_length))
    end
    pause(0.1)
end



function dydt = field(t,y)
%y(1)=P_x, y(2)=P_y
eqn1=y(2);    %v_x 
eqn2=-y(1);     %v_y
% eqn1=sin(y(2));    %v_x 
% eqn2=-y(1);     %v_y
dydt = [eqn1;eqn2];

end