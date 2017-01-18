%%================================================================================
%% This file is part of the bacterial envelope parameters
%%
%% (C)2016 Antonio Prestes Garcia <@>
%% For license terms see DESCRIPTION and/or LICENSE
%%
%% $Id$
##================================================================================

%{
   @title sim_growth

%}
function output = sim_growth(model, l0, ln, w, rho, h, s)
   switch model
      case 1
         output = sim_elongation1(l0, ln, w, rho, h, s);
         
      otherwise
	      fprintf('Invalid simulation function! Valid options are:\n');
	      fprintf('\t 1 for linear uptake\n');
	      fprintf('\t 2 for uptake proportional to surface area\n');
   end
   
   if exist('output', 'var'), sim_plot(output); end
end


%{
 	@title Simulates the bacterial cell envelope elongation and growth
%}

function output = sim_elongation1(l0, ln, w0, rho, steps, rodshape)
   %% ----- Handling missing parameters
   if ~exist('l0', 'var'), ln = 1.0; end
   if ~exist('ln', 'var'), ln = 2.0; end
   if ~exist('w0', 'var'), w0 = 0.5; end
   if ~exist('rho', 'var'), rho = 1.0; end
   if ~exist('steps', 'var'), steps = 100; end
   if ~exist('rodshape', 'var'), rodshape = true; end
   
   %% ----- Initialization
   rho = rho * 1000.0;              %% ----- SI, Kg/m3
   output = zeros(steps,7);         
   
   %% ----- Parameter validation
   if ~rodshape && ~(l0 == w0), fprintf("Spherical cells should have l0 == w0\n"); return; end
   
   w = w0;
   m0 = normrnd(rho, rho * 0.1) * b_volume(l0 * micrometer, w0 * micrometer);
   mn = normrnd(rho, rho * 0.1) * b_volume(ln * micrometer, w0 * micrometer);
   delta = (mn - m0) / steps;
      
   %% ----- Main simulation loop
   for i = 1:steps
      [s, v] = b_envelope(l0, w);
      m = m0 + normrnd(delta, delta * 0.1);    %% ---- Simple linear uptake with gaussian noise             
      l = b_length(m0, m, w0, w, l0);   
      output(i, :) = [ i, l0, v, s, (s/v), (m/femtogram), l];                  
      
      %% ----- Updating values for next step
      if ~rodshape, w = w0 = l; end
      m0 = m;         
      l0 = l;
   end
   %%BctPlotCellData(output); 
end

%{ 
   @title b.volume
   @description Calculates the bacterial envelope volume
   
   @param l The bacterial envelope length
   @param w The bacterial envelope width
   
   @return The volume

%}
function output = b_volume(l, w)
   output = ( (w^2 * pi/4) * (l-w) ) + (pi * w^3/6);
end

%{ 
   @title b.area
   @description Calculates the bacterial envelope area
   
   @param l The bacterial envelope length
   @param w The bacterial envelope width
   
   @return The volume

%}
function output = b_area(l, w)
   output = pi * w * (l - w) + pi * w^2;   
end

%{ 
   @title b.envelope
   @description Calculates the bacterial envelope area and volume
   
   @param l The bacterial envelope length
   @param w The bacterial envelope width
   
   @return The area and volume

%}
function [a, v] = b_envelope(l, w)
   a = b_area(l, w);
   v = b_volume(l, w); 
end

%{
   @title micrometer
   @description Returns the value of 1 micrometer in metters
   
   @return Returns the value of 1 micrometer in metters
%}
function v = micrometer()
	v = 1.0000e-006;
end

%{
 	Constant femtogram
 	Returns the value of 1 femtogram (10^-15 grams) in kilograms
    - "fg", femtogram, 10^-018 kilograms; 
%}
function v = femtogram
	v = 1.00e-018; 
end

%{ 
   @title b_length
   @description Calculate the updated length
 
   @param m0 The bacterial mass at time t
   @param m1 The bacterial mass at time t+1  
   @param w0 The bacterial width at time t
   @param w1 The bacterial width at time t+1
   @param l0 The bacterial lenght at time t
   
   @return The length at time t+1

%}
function output = b_length(m0, m1, w0, w1, l0)
   output = (w1^3 * m1 - w0^3 * m1 + 3 * w0^2 * l0 * m1) / (3 * m0 * w1^2); 
end

%{ 
   @title sim_plot
   @description Generate the plots for simulation output
 
   @param output The simulation output 
   
%}
function sim_plot(output)
   S = output;
	subplot(2,2,1); hold all;
   plot(S(:,2), S(:,3)); xlabel('Length (micrometers)'); ylabel('Volume (micrometers^3)'); title('Cell volume'); hold all; 
   subplot(2,2,2); hold all;
   plot(S(:,2), S(:,4)); xlabel('Length (micrometers)'); ylabel('Surface (micrometers^2)'); title('Cell surface'); hold all;
   subplot(2,2,3); hold all;
   plot(S(:,2), S(:,5)); xlabel('Length (micrometers)'); ylabel('Surface/Volume'); title('Cell S/V ratio'); hold all;
   subplot(2,2,4); hold all;
   plot(S(:,2), S(:,6)); xlabel('Length (micrometers)'); ylabel('Mass (femtograms)'); title('Cell mass');  
   %subplot(2,2,5); hold all;
   %plot(S(:,2), S(:,7)); xlabel('Length (micrometers)'); ylabel('Length (micrometers) estimated'); title('Cell Length');  
end
