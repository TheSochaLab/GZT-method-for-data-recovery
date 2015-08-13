% The Generalized Z-Transform (GZT) method is used to recover the 
% instantaneous gas exchange signal in flow-through respirometry systems.
% This is an improved version of the Z-transform (or Bartholomew) method.
% The code here enables a scientific user apply this method to their data.
% The details and rational for the method are published here: 
% 
% Pendar, H. and Socha, J.J. 2015. Estimation of instantaneous gas exchange in flow-
% through respirometry systems: A modern revision of Bartholomew's 
% Z-transform method, PLoS One. 
% 
% The method involves two steps: 
% 
% 1-	Calibration: 
% Before applying the method, the experimental setup must be 
% calibrated by finding the constant vector 'a' (Eq. 13 and 14). After 
% setting the flow rate, CO2 (or any other gas of interest) should be injected 
% into the respiratory chamber with a known pattern while recording output of the gas 
% analyzer. The input and output signals should be recorded in a text file 
% with three columns (arrays) of numerical data: time, concentration of the input gas, and 
% concentration of the output gas. This text file should be saved with name 
% of 'CalibrationData.txt' in the folder. After running the code and 
% choosing the calibration data file, the code will   
% find the calibration constants. The result will be saved in the file   
% 'Parameters_a.txt.
% 
% 2-	Instantaneous Gas Exchange Signal Recovery:
% After performing an experiment on the organism of interest, the raw gas exchange signal should be saved in a 
% text file with the name of ‘RawData.txt’ with two columns: time and gas concentration. After running this 
% code and choosing to recover the instantaneous signal, the code reads this 
% file and uses the calibration constants to recover the instantaneous gas 
% exchange signal. The result will be saved in a text file with the name of 
% 'Recovered_Data.txt’, which contains two columns of data: time
% and the recovered signal.
% 
% 
% 
% 
% 
%  H. Pendar August 2015
% 

%% Initializing:
close all
clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calibration
Calibrate = input('Would you like to calibrate (y/n)?      ','s');

if (Calibrate == 'y') || (Calibrate == 'Y') 
    % Loading files
    P=textread('CalibrationData.txt');  % reading the Calibration data
    time=P(:,1);  % time
    u_actual=P(:,2);  % actual input
    y=P(:,3);  % recorded signal in the output

    n=length(y);

    % Finding the parameters:

    % N could be any number less than the length of the impulse response
    N=230;  % number of parameters, should be chosen by trial and error

    % finding the parameters based on Eq. 14
    C=zeros(n-N+1,N);
    for i=1:N
        C(:,i)=y(i:n-N+i)';
    end

    a=inv(C'*C)*C'*u_actual(1:n-N+1);  % initial guess for the linear parameters
    figure
    plot(1:N, a)
    set(gca,'fontsize', 18)
    title('parameters (vector a in Eq.14)', 'FontSize',20)
    drawnow

    u_recovered=smooth(C*a,5,'moving'); % Input estimation based on Eq. 13. 'u_recovered' corresponds to variable 'u' in the paper.

    figure
    plot(time,u_actual)
    hold on
    plot(time(1:length(u_recovered)),u_recovered,'r');
    set(gca,'fontsize', 18)
    title('Blue: Actual input       Red: Estimation of the input', 'FontSize',20)
    drawnow

    save Parameters_a.txt a -ASCII   % Saves the obtained parameters in a text file

end

SignalRecovery = input('Would you like to recover the instantaneous signal (y/n)?      ','s');

if (SignalRecovery == 'y') || (SignalRecovery == 'Y')
    
    % Loading the files:
    a = textread('Parameters_a.txt');
    N = length(a);
    
    Data = textread('RawData.txt');   % The data should be saved as 'RawData.txt' in the same folder. It should contain two columns: time and the recorded gas exchange signal.
    time = Data(:,1);
    c = Data(:,2);
    n = length(c);
    
    C = zeros(n-N+1,N);
    for i=1:N
        C(:,i)=c(i:n-N+i)';
    end
    
    u_recovered = C * a;   % Instantaneous gas exchange calculation
    
    % noise filtering
    u_recovered = smooth(u_recovered,5, 'moving');  % A moving average is applied to smooth the results. By increasing the window size (which is arbitrarily 5 here), you will smooth data the data more, but the temporal accuracy will decrease. 
    
    Threshold = 0.95;  % The threshold can be found by running this code over a pure noise signal of the system. Any recovered signal below the threshold will be significanly smoothed out.
    u_recovered(u_recovered < Threshold) = smooth(u_recovered(u_recovered < Threshold),20, 'moving');
    
    figure
    
    plot(time, c, 'linewidth', 2)
    
    hold on
    m=min(length(time), length(u_recovered));
    plot(time(1:m), u_recovered(1:m), 'r', 'LineWidth', 2)
    
    xlabel('Time', 'FontSize', 20)
    ylabel('Gas Concentration', 'FontSize', 20)
    
    title('Blue: Raw recorded signal      Red: Instantaneous recovered signal', 'FontSize',20);
    set(gca,'fontsize', 18)
    
    % Save the recovered instantaneous signal in a text file:
    Output = [time(1:m) , u_recovered(1:m)];
    save Recovered_Data.txt Output  -ASCII
    
end
