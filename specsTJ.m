% 1D SPECS code readapted from the Yuhai Tu's 2D version, originally by
% EDS.
% Changing this 
% Test with this command: 
% specs1D(0,200,150,100,100,0.005,0);

function specsTJ(ligandleft,ligandright,preadaptC,total_t,num_cell,kr,flag)
third_argument=preadaptC
display_progress=1;    %set to 0 not to printout anything, 1 if want printout

% flag is 0=e.coli,1=salmonella; 2=e.coli with inst. tumbling
% e.g.  (0,200,150,100,100,0.005,0);   for preadapted to 150 and total time 100
% IMPORTANT:  default in SPECS for e.coli is: kr = 0.005
% or  specs_expo_linear(0,0,2000,100,100,0.05,1)  
%       for char lengthscale 2000 in expo ligand, kr=0.05, 100 salmonellas

% simulation time total_t is in seconds
% third_argument is either "characteristic lengthscale of the gradient in um"
%  (for exponential gradients) or (linear case) just the initial ligand

% to show histograms make =1, else=0
showdistributions=1;
interval_between_frames=0.001;  % e.g. 0.1, 0.01
% is above=1 then a "movie" is shown, with that interval between frames
% seems not to plot intermediates unless I include this?

showmeans=1;
% if==1, end by plotting the history of means, so no need to call plot prog

uniforminitial=1;
% set to 1 for uniform initial distribution, zero for not

middleinitial=0;
% set to 1 for delta initial distro, zero for not
% centered in the middle unless exponential ligand (in which case, computed)

linear_ligand=1;
% set to 1 for linear ligand, zero for not

exponential_ligand=0;
% set to 1 for exponential ligand, zero for not

% Sort out the conflict between total_x and channellength
total_x=1;
        
if display_progress==1
     if flag==0
        display('e.coli');
     end
     if flag==1
        display('salmonella');
     end
     if flag==2
        display('e.coli with instantaneous tumbling');
     end
     if linear_ligand==1
        display('linear ligand')
     end
     if exponential_ligand==1
        display('exponential ligand')
     end
     if uniforminitial==1
        display('uniform initial condition')
     end
     if middleinitial==1
        display('delta initial condition')
     end
end 

length_channel = 1000; % unit: um; default, not used for expo gradient!
                   %because right now I want to analyze only initial behavior
dx = 0.01;                % must make tiny, else too much round-off

displaymethylation=0;
% set to zero for not displaying methylation levels at end

dt = 0.001;             %unit: s. time step for Euler solving; MUST be small

save_interval=1;
% s; every how many seconds to take a snapshot and save data
sample_dt = save_interval/dt;
% every how many samples, in dt units, to take a snapshot and save data
times=[];     % to store vector of times at which took snapshots, for plotting

% create a filename for the saved data, for expo ligand:
if exponential_ligand==1
 rootname = 'specs_saved';
 underscore='_';
 extension = '.mat';
 filename = [rootname, underscore, num2str(kr), underscore, num2str(third_argument),extension];
end

% create a filename for the saved data, for linear ligand:
if linear_ligand==1
 rootname = 'specs_saved';
 underscore='_';
 extension = '.mat';
 filename = [rootname, underscore, num2str(ligandleft), underscore, num2str(ligandright),extension];
end;

number_steps = total_t/dt;    % # time points at which Euler will be computed

% length_channel is the length of the channel
% dx is the width of the discrete grid on which cells are placed
% so that of length=200 and dx=10, there are gridpoints at positions 1...20 
%     I should change to "continuous" position and computing ligand as needed
%     but of course that means more computing rather than a simple table lookup

% should also define a bigger "DX" for histogram, to filter out noise 

allmeans=[];  % to collect mean values (eds)

if flag==0|2    % specific parameters for e.coli:
 KI_TAR_MEASP = 18.2;
	% um; dissociation constant of MeAsp to inactive Tar receptor
 KA_TAR_MEASP = 3000;
	% um; dissociation constant of MeAsp to active Tar receptor
 N_TAR = 6;
	% number of Tar receptors dimers in a complex
 ALFA_TAR_MEASP = 1.7;
	% parameter in equation (2) of Jiang et al
 M0_TAR_MEASP = 1;
	% parameter in equation (3)
% KR = 0.005;
 KR = kr;
	% 1/s;  linear rate for methylation process
 KB = KR;
	% 1/s;  linear rate for demethylation process
 H = 10.3;
	% Hill coefficient of CheY-P response curve
 A_0 = 0.5;
	% steady state activity value
 RUN_VELOCITY = 16.5;
        %um/s;  average run velocity
end % specific to e.coli

if flag==1   % specific parameters for salmonella:
 KI_TAR_MEASP = 27;
	% um; dissociation constant of MeAsp to inactive Tar receptor
 KA_TAR_MEASP = 27000;
	% um; dissociation constant of MeAsp to active Tar receptor
 N_TAR = 2;
	% number of Tar receptors dimers in a complex
 ALFA_TAR_MEASP = 1.7;
	% parameter in equation (2) of Jiang et al
 M0_TAR_MEASP = 1;
	% parameter in equation (3)
% KR = 0.005;
 KR = kr;
	% 1/s;  linear rate for methylation process
 KB = KR;
	% 1/s;  linear rate for demethylation process
 H = 10.3;
	% Hill coefficient of CheY-P response curve
 A_0 = 0.5;
	% steady state activity value
 RUN_VELOCITY = 40.0;
        %um/s;  average run velocity
end % specific to salmonella

% common parameters:
RUN_TIME = 0.8;
	%s;   average run time
TUMBLE_TIME = 0.2;
	%s;   tumble time (deterministic)
TUMBLE_STEP_MAX = TUMBLE_TIME / dt;
        % for how many [Euler] steps we tumble (deterministic)

if display_progress==1
     if flag==2
        TUMBLE_TIME = 0.0;
     end
end

% linear ligand concentration:
if linear_ligand==1
  number_mesh_points = length_channel/dx;
% initialize ligand concentration at zero (will make gradient below)
  ligand_1d = zeros(1,number_mesh_points);
  for i = 1 : 1 : number_mesh_points
   ligand_1d(i) = (ligandright-ligandleft)/number_mesh_points * i + ligandleft;
  end
  initial_ligand = third_argument;
end

% could also write as "else":
if exponential_ligand==1
  maxExpectedVd=20;
  length_channel = maxExpectedVd*total_t; % as per Tom's email
  initial_ligand = sqrt(KI_TAR_MEASP*KA_TAR_MEASP);
  length_scale=third_argument;
  number_mesh_points = length_channel/dx;
% initialize ligand concentration at zero (will make gradient below)
  ligand_1d = zeros(1,number_mesh_points);
  C=3500;  % um^2/s
  G_C=sqrt(2*ALFA_TAR_MEASP*KR*A_0/(H*C));
  v_d = C*G_C / (1 + length_scale *G_C );
  initial_position = length_channel/2 - (total_t/2)*v_d;
%check with F that (...) in deno
  for i = 1 : 1 : number_mesh_points
   ligand_1d(i) =initial_ligand*exp( (i*dx - length_channel /2) / length_scale);
  end
end

%  make sure to use a very thin mesh, else there is too much round-off,
%    since anyone between two mesh points is assigned value at right!
%    -- otherwise, evaluate function, instead of table lookup

% i.e. use ligand_1d, which is a vector, instead of repeating it at each time, 
%   because we are assuming constant.
%   but: if I wanted to have it varying, then use a structure like this:
%      ligand_all_times = repmat(ligand_1d, [1 1 number_steps]);
%   this makes a 3-d matrix in which for each t in 1:1:number_steps, we have
%      the 1d ligand; of course in this example, it is still constant

%% INITIALIZATION

x = zeros(num_cell, 1);
% means: x(i) = position of cell i

% this formula gives the steady-state methylation level corresponding to the initial_ligand

L=initial_ligand;
initial_methylation = M0_TAR_MEASP - (1/ALFA_TAR_MEASP)*(1/N_TAR*log(1/A_0-1)+(log((1+ L/KA_TAR_MEASP)/(1+L/KI_TAR_MEASP))));

%             % debug test if get 0.5
%             m = initial_methylation
%             a = 1 ./ (1 + exp(N_TAR .* ((ALFA_TAR_MEASP .* (M0_TAR_MEASP - m)) - (log((1 + L./KA_TAR_MEASP)./(1 + L./KI_TAR_MEASP))))))


m = initial_methylation*(zeros(num_cell, 1) + 1);
% m(i) = methylation level of cell i, all equal "initial_methylation" at t=0

% this is a memory of what the cell was doing before we loop on update
tumble_prev = zeros(num_cell, 1);
% tumble_prev(i) takes values 0=tumble or 1=run; start with "tumbling" at t=0

tumble_counter = zeros(num_cell, 1) + 1;
% tumble_counter(i) = how long cell i has been tumbling; starts at "1" at t=0

direction = zeros(num_cell, 1);                 %Cell's run direction, 0-360;
% direction(i) = run direction of cell i in degrees; will be -1 or 1 in 1d

L = zeros(num_cell, 1);
% L(i) = ligand concentration at the location of each cell i

%num_save = total_t / dt_save;
num_save = total_t / save_interval;
% how many snapshots to save

% making sparse matrix, in large channel most won't be used in small time
save_x = sparse(num_cell, num_save);
%save_x = zeros(num_cell, num_save);
% save_x(i,t) = where cell i was at time t

cell_distributions_all = sparse(number_mesh_points, num_save);
%cell_distributions_all = zeros(number_mesh_points, num_save);
% cell_distribution_time(k,t) = how many cells in bin "k" at time t

save_m = sparse(num_cell, num_save);
%save_m = zeros(num_cell, num_save);
% save_m(i,t) = methylation level of cell i at sample time t

save_a = sparse(num_cell, num_save);
%save_a = zeros(num_cell, num_save);
% save_a(i,t) = activity level of cell i at sample time t

%initialize the location of each cell
%  -- if nothing said, all x(i) are = 0

% code for case when they start uniformly distributed:
if uniforminitial==1
    for i = 1 : 1 : num_cell    
    %Uniform distribution:
     temp = round (i * (length_channel * total_x) / (num_cell + 1));
      x(i) = ceil(temp/total_x);
    end
end

% code for case when they start all in a delta function:

midpositionx=round((1/2)*length_channel);  % default

if exponential_ligand
   midpositionx = initial_position;
end

if middleinitial==1
    for i =1:1:num_cell
      x(i) = midpositionx;
    end
end
        
if showdistributions==1
          figure(1)
end
% just to make sure I plot in window #1

%% RUN SIMULATION NOW

% start main loop
for t = 1 : 1 : number_steps
                  %debug
                  %if  mod(t, 10000)==0
                  %  fprintf('step=%d out of %d\n',t,number_steps)
                  %end
    s = ligand_1d;
%                          %debug
%                          x
    for n = 1 : 1 : num_cell
        temp1 = ceil(x(n)/dx);
% temp1 has the position of cell n at the current time, rounded up to grid point
% I am not sure how x(n) could ever have been zero, but if it was, rounding-up
%    would not work (index=0 out of bounds), so mapping '0' to '1' grid point:
        if temp1 == 0
            temp1 = 1;
        end % correction for zero grid point
% measure ligand concentration at current position (grid point):
       L(n) = s(temp1);
    end % of computing ligands for each cell; L is a vector now

% compute activity  (equations 1 and 2 in Jiang et al)
% remember that m(i) = methylation level of cell i, and these are vector ops
% the vector "m" is updated, and a vector a of activities is created

    a = 1 ./ (1 + exp(N_TAR .* ((ALFA_TAR_MEASP .* (M0_TAR_MEASP - m))...
        - (log((1 + L./KA_TAR_MEASP)./(1 + L./KI_TAR_MEASP))))));

% compute methylation (equation 3 in Jiang et al)

    m = m + (KR .* (1 - a) - KB .* a) .* dt;

% probability of change from run to tumble in interval of length dt:
%                                (equations 4 and 5 in Jiang et al)
% note that this is a vector, listing the probability for each cell

  %  p = (dt / RUN_TIME) .* L./ligandright;%ones(num_cell,1);%(a./A_0).^H;
    p =  (1-(sin(2*pi*2.5.*x/1000)+1)/2).^2;
                  %debug :print probabilities
                  %if p >0
                  %   if p < 1
                  %      display('****************')
                  %      p
                  %   end
                  %end
% loop to decide if to tumble or run, and then update position

    for i = 1 : 1 : num_cell

% for each cell,
% if was tumbling already, and if the tumbling time is < TUMBLE_STEP_MAX
        if tumble_prev(i) == 0
              if tumble_counter(i) < TUMBLE_STEP_MAX
% then stay tumbling: update tumble memory and counter of tumbling time:
                  tumble_prev(i) = 0;
                  tumble_counter(i) = tumble_counter(i) + 1; 
% but if we got to the end of the tumble period,
              else
% then change mode to run:
                  tumble_prev(i) = 1;

% and make direction = 1 or -1 with prob 1/2:
                  direction(i) = 2*ceil(rand(1)-0.5)-1;
%                         %debug
%                         display('changed direction to:')
%                         direction(i)
              end     % of "if not at end of tumbling, else"

% if, instead, the cell was running,
        else
% then with probability p(cell), start tumbling:
            r = rand(1);    %Get a random number between 0 and 1
            if r < p(i)
%                          %debug
%                          p
%                          display('start tumbling')
                tumble_prev(i) = 0;     % tumble
                tumble_counter(i) = 1;  % start again the tumble step counter	
% and with probability 1-p continue running:
            else
                tumble_prev(i) = 1;    %run
            end    % of "generating probability of tumble"
        end        % of "if tumble prev=0" i.e. "if was tumbling"

% still for this cell, if I will keep running (note "tumble_prev" has been updated)
        if tumble_prev(i) == 1
% then update the position
           x(i) = x(i) +  direction(i)*RUN_VELOCITY * dt;
        end   % of updating position if running

% next, detect the boundary and deal with the boundary effect.
% for 1d, just make a reflection, as in theoretical Neumann condition 

       if x(i) < 0
           x(i) = -x(i);    % reflect
           direction(i) = -direction(i); 
%                          %debug
%                          display('***** hit boundary:')
%                          x(i)
%                          direction(i)
      elseif x(i) > length_channel
           x(i) = 2*length_channel - x(i);  % = total_x-(x(i)-total_x), i.e. reflection
           direction(i) = -direction(i);
%                          %debug
%                          display('***** hit boundary:')
%                          x(i)
%                          direction(i)
       end % of boundary reflection


    end       % of main loop for tumbling or running and updating position

%  now collect and save data, provided
%  current time t*dt is a multiple of sampling time:

    if  mod(t, sample_dt)==0
        if display_progress==1
             % display dots while thinking...:
             fprintf('.');
        end
        save_x(:, t/sample_dt) = x;
        save_m(:, t/sample_dt) = m;
        save_a(:, t/sample_dt) = a;
        times = [times,t/sample_dt];

% (saved vectors x,m,a (one entry per cell) at the corresponding sample time)

% compute in which box does the coordinate x lie [rounding up so never 0]:

        cell_distribution_current = zeros(number_mesh_points,1);
        for i= 1 : 1 : num_cell
            temp1 = ceil(x(i)/dx);
                if temp1 == 0
                   temp1 = 1;
                end
%                      and now simply add +1 to the corresponding pixel:
            cell_distribution_current(temp1)=cell_distribution_current(temp1)+1;
        end      % of loop that updated current cell distribution vector

% store the current distro in the array cell_distributions_all:

    cell_distributions_all(: , t/sample_dt) = cell_distribution_current;

% storing them, but not yet doing anything with the result!

% plot them, with pause to be able to visualize
      if showdistributions==1
          plot(cell_distribution_current)
          pause(interval_between_frames);
      end

[howmany,dummy]=size(cell_distribution_current);
range=1:howmany;
meanvalue=sum(range*cell_distribution_current);
% need to correct units, though:
meanvalue = meanvalue*dx/num_cell;
% now collect
allmeans=[allmeans,meanvalue];

    end   % of saving at sample times
                  %debug : print values of x (cell positions) at all times
                  %if p >0
                  %   if p < 1
                  %      display('****************')
                  %      x
                  %   end
                  %end
end  %  of main loop with data collection at sample times 

if display_progress==1
   fprintf('\n')
   display(['saving data in file: ',filename]);
end

%                        %debug
%                        x

%print vector of means:
%allmeans

% the data at the last iteration is in:  cell_distribution_current
% saving the structure with all the frames as well as methylation, etc
%
% uncomment if want to save data:
%
%savefile = 'specs_results_saved.mat';
save(filename,'save_x','save_m','save_a','times','allmeans','cell_distributions_all','ligand_1d');

if displaymethylation==1
% show methylation levels at end
   format long
   save_m(:,end)'
   format short
end

if showmeans==1
   figure(3)
   plot(times,allmeans)
end

% uncomment if want to print:
%
%allmeans

