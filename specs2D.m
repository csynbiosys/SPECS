% changing time to very long to be able to get a steady-state

function [save_x save_y save_m save_a F]=specs2D(left,right,initial_methylation)

fprintf('.')

% (eds) added "marginal=" so that can pass out the final marginal distribution
%  also pass arguments for left and right ligand concentrations
%  also added the initial methylation level as input; was '1' in original code
%  examples of initial levels:
  %  2.081274522478973   for 100 constant
  %  2.423207034223378   for 200 constant
  %  2.535452158288703   for 250 constant
  %  2.704033304437643   for 350 constant
  %  2.738302210901962   for 375 constant
  %  3.1981   for 1000 constant
  %  3.4692       2000
%figure(1)

%Rewrote layout so that size 20 units includes 0 & 20, ie 20 spaces, 21
%nodes.


%This a 2D SPECS program.
% (eds) adding, '0' not so show movie,'1' to show
showmovie=1;

% (eds) adding, to show histograms:
showdistributions=0;

% (eds): ligand at boundaries:
ligandleft=left;
ligandright=right;
%hold on

uniforminitial=1;
% set to zero for uniform initial distribution, zero for not

middleinitial=0;
% set to zero for delta at middle initial distribution, zero for not

displaymethylation=0;
% set to zero for not displaying methylation levels at end


%% Set parameters
%Simulation information 
%total_t = 1600;         %unit: s. The total simulation time
total_t = 2400;         %unit: s. The total simulation time
dt = 0.1;               %unit: s. Time step
%dt_save = 20;           %unit: s. The period to save data
dt_save = 20;           %unit: s. The period to save data

%total_x = 800;          %unit: um. Channel length
%dx = 5;                 %unit: um
%total_y = 400;          %unit: um. Channel width
%dy = 5;                 %unit: um

% above are original values
total_x = 1000;          %unit: um. Channel length
dx = 10;                 %unit: um

total_y = 1000;          %unit: um. Channel width
dy = 10;                 %unit: um

%num_cell = 2000;     	%Number of the cells
num_cell = 10000;     	%Number of the cells

allmeans=[];  % to collect mean values (eds)

%Define Position Grid Matrix for interpolation, etc
[Xgrid,Ygrid]=meshgrid([0:dx:total_x],[0:dy:total_y]);
% Xgrid=Xgrid';
% Ygrid=Ygrid';


%Pathway related parameters
KI_TAR_MEASP = 18.2;		%unit: um. Dissociation constant of MeAsp to inactive Tar receptor
KA_TAR_MEASP = 3000;		%unit: um. Dissociation constant of MeAsp to active Tar receptor
N_TAR = 6;					%Number of Tar receptors dimers in the complex

ALFA_TAR_MEASP = 1.7;		%A parmeter in equation (2)
M0_TAR_MEASP = 1;			%A parmeter in equation (3)

KR = 0.005;					%unit: 1/s. Linear rate for methylation process
KB = KR;					%unit: 1/s. Linear rate for demethylaion process

H = 10.3;					%Hill coefficient of CheY-P response curve
A_0 = 0.5;					%Activity value in steady state

%Run & Tumble related paramerers
RUN_VELOCITY = 16.5;		%unit: um/s. Average run velocity
RUN_TIME = 0.8;				%unit: s. Average run time
TUMBLE_TIME = 0.2;			%unit: s. Constant tumble time 
TUMBLE_STEP_MAX = TUMBLE_TIME / dt;
ROT_DIF_CONST = 30;			%Average directional change in 1s

%% Get the Ligand Concentration Profile
%Generate the ligand concentration profile or load the file contain the
%ligand concentration information (make sure the dimensions of the matrix 
%are consistant with the parameters defined above).

%Method 1, generate a ligand concentration profile:
ligand_2d = zeros(total_x/dx+1, total_y/dy+1);
for i = 1 : 1 : total_x/dx+1
    for j = 1 : 1 : total_y/dy+1
%       ligand_2d(i, j) = 200 * (i/(total_x/dx)) + 400;
% make it steeper:
        ligand_2d(i, j) = (ligandright-ligandleft) * (i/(total_x/dx+1)) + ligandleft;
    end
end
% (eds comment:) repmat is like a Krokecker product with identity; I think
% that this is telling it to make one of these for each time instant, i.e.
% a constant ligand profile
ligand = repmat(ligand_2d, [1 1 total_t/dt]);

% %Method 2, define your ligand concentration and save in a file. Load it here:
% % load data_ligand.mat

%% Declaration and Initialization
x = zeros(num_cell, 1);                         %The position of each cell
y = zeros(num_cell, 1);
m = initial_methylation*(zeros(num_cell, 1) + 1);                     %The methylation level of each file
mold = zeros(num_cell, 1);                      %0 = tumble; 1 = run.
tumble_counter = zeros(num_cell, 1) + 1;        %The duration of a tumble
direction = zeros(num_cell, 1);                 %Cell's run direction, 0-360;
L = zeros(num_cell, 1);                         %Ligand concentration at the location of each cell

num_save = total_t / dt_save;                   %The total steps of saving
save_x = zeros(num_cell, num_save);             %Save the cell's positions
save_y = zeros(num_cell, num_save);
save_cell_dis = zeros(total_x/dx+1, total_y/dy+1, num_save); %Save the cell distribution
save_m = zeros(num_cell, num_save);            %Save the methylationlevey level
save_a = zeros(num_cell, num_save);            %Save the activity of each cell

%initialize the location of each cell
% (eds comment:) if nothing said, all x(i) and y(i) are = 0

%V - This needs to be revisited to check it works with revisions
if uniforminitial==1
    for i = 1 : 1 : num_cell    
    %Uniform distribution:
     temp = round (i * (total_x * total_y) / (num_cell + 1));
      x(i) = ceil(temp/total_y);
      y(i) = mod(temp, total_y);
    end
end

%V - Revised to remove pointless loop
if middleinitial==1
    midpositionx=round((1/2)*total_x);
    x(:)=midpositionx;
    midpositiony=round((1/2)*total_y);
    y(:)=midpositiony;
end
        
%% SPECS

%Main Time Loop
jj=1;
for t = 1 : 1 : total_t/dt
    s = ligand(:, :, t);    %Get the ligand profile at current time point

    %V - This appears to check whether the bacteria locations are at the lower
    %boundary..... but they are initialized to that boundary anyway?
    
    %This also snaps the cell position to a grid for the ligand
    %concentration.... it should be interpolated to the true location.
    %changed to the following:
    L=interp2(Xgrid,Ygrid,s',x,y);
    
    %Get the kinase activity (equation 1 and 2)
    a = 1 ./ (1 + exp(N_TAR .* ((ALFA_TAR_MEASP .* (M0_TAR_MEASP - m))...
        - (log((1 + L./KA_TAR_MEASP)./(1 + L./KI_TAR_MEASP))))));
    %Update the mathylation level (equation 3)
    m = m + (KR .* (1 - a) - KB .* a) .* dt;
    %Get the probility density of CW (equation 4 and 5)
    p = (dt / RUN_TIME) .* (a./A_0).^H;
    
%   Rework this section to not loop as much as possible  

    %Find New Tumblers
    f=find(mold);
    g=find(rand(size(f))<p(f));
    mold(f(g))=0;
    tumble_counter(f(g))=0; %initialize tumble count to zero since it will be updated in this same period
    %Run Update
    f=find(mold);
    direction(f)=mod(direction(f) + randn(size(f)) * ROT_DIF_CONST * dt.^0.5,360); %Rotational Diffusion update
    x(f) = x(f) + RUN_VELOCITY * dt * cos(direction(f) * pi / 180); %Update the position
    y(f) = y(f) + RUN_VELOCITY * dt * sin(direction(f) * pi / 180);
    g=find(x(f)<0); %check boundary enforcement
    x(f(g))=0;
%     direction(f(g))=floor(direction(f(g)))/180* 180 + 90;  %V - do they change direction??
    g=find(x(f)>total_x); %check boundary enforcement
    x(f(g))=total_x;
%     direction(f(g))=floor(direction(f(g)))/180* 180 + 90;  %V - do they change direction??
    g=find(y(f)<0); %check boundary enforcement
    y(f(g))=0;
%     direction(f(g))=(1 - floor(direction(f(g))/270))* 180;  %V - do they change direction??
    g=find(y(f)>total_y); %check boundary enforcement
    y(f(g))=total_y;
%     direction(f(g))=floor(direction(f(g))/90)* 180;  %V - do they change direction??
    %Tumble Update
    f=find(and(~mold,tumble_counter < TUMBLE_STEP_MAX)); %find tumbling cells with time still to tumble
    tumble_counter(f) = tumble_counter(f) + 1; 
    f=find(and(~mold,tumble_counter >= TUMBLE_STEP_MAX)); %find tumbling cells now done tumbling
    mold(f)=1;
    direction(f)=rand(size(f))*360;

% figure(3)
% plot(x,y,'.')
% %     axis([0,total_x,0,total_y])
% pause

    %Save data
    if  mod(t*dt, dt_save)==0
        % (eds) i.e., if current time t*dt is a multiple of sampling time:
        save_x(:, t*dt / dt_save) = x;
        save_y(:, t*dt / dt_save) = y;
        save_m(:, t*dt / dt_save) = m;
        save_a(:, t*dt / dt_save) = a;
        
        cell_dis=hist3([x,y],{[0:dx:total_x]',[0:dy:total_y]'});
        save_cell_dis(:, :, t*dt / dt_save) = cell_dis;
% (eds) above line adds to the structure save_cell_dis the current distro
%       -- i.e. the first two coordinates are the distros and the last is
%          the time frame

%Check the real time cell distribution 
% (eds) added conditional showing pictures:
        if showmovie==1
             figure(1)
             ftemp = round( (cell_dis-min(min(cell_dis)))./(max(max(cell_dis)) - min(min(cell_dis))) *64 ) ;
             ftemp_show = imresize(ftemp, [total_x/2 total_y/2]);
             imshow(ftemp_show, jet,'InitialMagnification',55);
             imagesc(ftemp_show)
             colormap(jet)
             axis equal
             title(['t=', num2str(t*dt), 's']);
             colorbar;
             axis tight
             pause(0.1);
             F(1,jj)=getframe;
             jj=jj+1;

        end
% (eds adding command to add along length:)
        marginal=sum(cell_dis');
        if showdistributions==1
          figure(2)
          plot(marginal)
          pause(0.01);
        end
        howmany2=size(marginal);
        howmany=howmany2(1,2);
        range=1:howmany;
        meanvalue=sum(range.*marginal);
        allmeans=[allmeans,meanvalue];
% (eds commands above) prints average displacement
    end
% (eds) ; now end data collection at sample times 
end
% here is the data at the last iteration:  cell_dis
savefile = 'data.mat';
% (eds) saving the structure with all the frames as well as methylation, etc
%
% uncomment if want to save data:
%
%save(savefile, 'save_x', 'save_y', 'save_m', 'save_a', 'save_cell_dis');
if displaymethylation==1
% show methylation levels at end
   format long
   save_m(:,end)'
   format short
end
% uncomment if want to print:
%
%allmeans
%
