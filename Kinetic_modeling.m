% Simutaneously fit data by Casper with double-conformation model 
% using non-linear least-squares fit
% Fit two data sets simutaneously
clc; clearvars; close all
% Set figure renderer
set(0, 'DefaultFigureRenderer', 'painters');

% Input k1, k2 to fit k3~k6
% k1 and k2 can both be [] to fit all 6 parameters without restraints
% I found that restraining k1 and k2 results in less error in the fit
k1 = 102; k2 = 1.1;
% Perform bootstrapping? (Y/N)
btrp = 'N';
% Automatically save figures? (Y/N)
save = 'N';

% Load data. (C8/ S6/ S8/ Sr/ Ca)
% Here I load the data of Ca2+-PSII at pH 8
[Ds, Dd, RNs_fit, RNd_fit] = loadData('C8');

% Concatenate Ys and Yd functions and data 
% to simutaneously fit the two data sets
fun = cctnt_fun(k1, k2);
[xdata, ydata] = cctnt_data(Ds, Dd);

%%% Non-linear curve fit
% First get x0 and lower bound. The number of inputs 
% depends on if k1 and k2 are empty
[x0, lb, ub] = getx0lb(k1, k2);
% Fit data and get residuals compared to simple exponential fit
x = lsqcurvefit(fun, x0, xdata, ydata, lb, ub);
[resnorms, resnormd] = res(Ds, Dd, x, k1, k2);

% Perform bootstrapping for error estimation
switch btrp
    case 'Y' 
        fit = @(x,y) lsqcurvefit(fun, x0, x, y, lb, ub);
        [~,bootstat] = bootci(5000,{fit, xdata, ydata});
end

% Print coefficients and residuals
format short g; 
if isempty(k1) && isempty(k2)
    disp('Fitted coefficients (k1~k6): '); disp([x(5:6) x(1:4)])
else
    disp('Fitted coefficients (k3~k6): '); disp(x)
end
disp('ResNorm of Ys = '); disp(resnorms-RNs_fit)
disp('ResNorm of Yd = '); disp(resnormd-RNd_fit)

% Plot data, fitted curves, and save figures
t = [linspace(0, 1.5, 1000) linspace(1.5, 12, 1000)];
[Ys, Yd] = getY(t, x, k1, k2);
myplot(t,Ys,Ds(:,1),Ds(:,2),'Single', save)
myplot(t,Yd,Dd(:,1),Dd(:,2),'Double', save)

% Load data
function [Ds, Dd, RNs_fit, RNd_fit] = loadData(data)
    switch data
    case 'C8'
        Ds = load('Ca_8.6_single.csv'); Dd = load('Ca_8.6_double.csv');
        RNs_fit = 0.0427; RNd_fit = 0.0377;
    case 'S6'
        Ds = load('Sr_6.0_single.csv'); Dd = load('Sr_6.0_double.csv');
        RNs_fit = 0.0651; RNd_fit = 0.0675;
    case 'S8'
        Ds = load('Sr_8.3_single.csv'); Dd = load('Sr_8.3_double.csv');
        RNs_fit = 0.0682; RNd_fit = 0.0884;
    case 'Sr'
        Ds = load('Sr_S2_single.csv'); Dd = load('Sr_S2_double.csv');
        RNs_fit = 0.0709; RNd_fit = 0.1339;
    case 'Ca'
        Ds = load('Ca_S2_single.csv'); Dd = load('Ca_S2_double.csv');
        RNs_fit = 0.1740; RNd_fit = 0.0327;
    end
end
% Concatenate functions
function fun = cctnt_fun(k1, k2)
    if isempty(k1) && isempty(k2)
        funs = @(x, xdata) M2(xdata, 3.9, x, 'Single');
        fund = @(x, xdata) M2(xdata, 3.9, x, 'Double');
    else
        funs = @(x, xdata) M2(xdata, 3.9, [x k1 k2], 'Single');
        fund = @(x, xdata) M2(xdata, 3.9, [x k1 k2], 'Double');
    end
    fun = @(x, xdata) [funs(x, xdata(:,1)), fund(x, xdata(:,2))];
end
% Concatenate data
function [xdata, ydata] = cctnt_data(Ds, Dd)
    if length(Ds) >= length(Dd)
        diff = length(Ds) - length(Dd);
        Dd = [Dd; zeros(diff, 2)]; 
    elseif length(Dd) > length(Ds)  
        diff = length(Dd) - length(Ds);
        Ds = [Ds; zeros(diff, 2)];
    end
    xdata = [Ds(:,1) Dd(:,1)];
    ydata = [Ds(:,2) Dd(:,2)];
end
% Get xo and lower/upper bounds depending on number of fitting parameters
function [x0, lb, ub] = getx0lb(k1, k2)
    if isempty(k1) && isempty(k2)
        x0 = [87.7 35.6 1    1    102  1.1 ];
        lb = [0    0    0    0    0    0   ];
        ub = [5000 500  500  500  inf  inf ];
    else
        x0 = [87.7 35.6 1    1   ];
        lb = [0    0    0    0   ];
        ub = [5000 500  500  5000];
    end
end
% Get Y values for fitted curves
function [Ys, Yd] = getY(t, x, k1, k2)
    if isempty(k1) && isempty(k2)
        [Ys, Yd] = M2(t, 6.4, x, 'Both');
    else
        [Ys, Yd] = M2(t, 6.4, [x k1 k2], 'Both');
    end
end
% Calculate residuals compared to simple exponential fittings
function [resnorms, resnormd] = res(Ds, Dd, x, k1, k2)
    if isempty(k1) && isempty(k2)
        Ys = M2(Ds(:,1), 6.4, x, 'Single');
        Yd = M2(Dd(:,1), 6.4, x, 'Double');
    else
        Ys = M2(Ds(:,1), 6.4, [x k1 k2], 'Single');
        Yd = M2(Dd(:,1), 6.4, [x k1 k2], 'Double');
    end      
    resnorms = sum((Ds(:,2) - Ys).^2);
    resnormd = sum((Dd(:,2) - Yd).^2);
end
% Plot data, fitted curves, and save figures
function myplot(x1, y1, x2, y2, mode, save)
    width = 3;     % Width in inches
    height = 2.5;  % Height in inches
    alw = 0.75;    % AxesLineWidth
    fsz = 11;      % Fontsize
    lw = 1;        % LineWidth
    msz = 2.3;     % MarkerSize
    
    figure()
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
    set(gca, 'FontSize', fsz, 'LineWidth', alw); %<- Set properties
    plot(x2,y2,'o','color',[0,0,0]+0.5,'MarkerFaceColor',[0,0,0]+0.5,...
        'MarkerSize', msz)
    hold on
    plot(x1,y1,'LineWidth',lw,'color','black')
    % myplotStyle()
    xlabel('H_2^{ 18}O Incubation Time (s)')
    ylabel('Normalized O_2 yield')

    ylim([0 1.1])
    xlim([-0.1 3.1])
    
    % Set Tick Marks
    set(gca,'XTick',0:0.5:3);
    set(gca,'YTick',0:0.2:1.2);
    set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
    set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
    % set(gca,'TickDir','out');
    
    % Save plot
    % Here we preserve the size of the image when we save it.
    switch save
        case 'Y'
            set(gcf,'InvertHardcopy','on');
            set(gcf,'PaperUnits', 'inches');
            papersize = get(gcf, 'PaperSize');
            left = (papersize(1)- width)/2;
            bottom = (papersize(2)- height)/2;
            myfiguresize = [left, bottom, width, height];
            set(gcf,'PaperPosition', myfiguresize);
            switch mode
                case 'Single'
                    print('C8 Single.png','-dpng','-r300');
                case 'Double'
                    print('C8 Double.png','-dpng','-r300');
            end
    end
end
% Fit data using the double-conformation model
function [Ys, Yd] = M2(t, r, k, mode)
    % Unpack k to k3, k4, k5, k6, k1, k2
    k = num2cell(k);
    [k3, k4, k5, k6, k1, k2] = k{:};

    % Eigenvalues of k matrix
    k56 = k5+k6;
    kT1 = k1+k2+k3+k4+k56;
    kT2 = k1+k3+k56;
    kT3 = k2+k4+k56;

    DL1 = sqrt(kT1^2-4*((k1+k2)*(k3+k4)+(k1+k2)*k5+(k3+k4)*k6));
    DL2 = sqrt(kT2^2-4*(k1*k3+k1*k5+k3*k6));
    DL3 = sqrt(kT3^2-4*(k2*k4+k2*k5+k4*k6));

    L1p = 0.5*(-kT1+DL1); L1m = 0.5*(-kT1-DL1);
    L2p = 0.5*(-kT2+DL2); L2m = 0.5*(-kT2-DL2);
    L3p = 0.5*(-kT3+DL3); L3m = 0.5*(-kT3-DL3);

    % Analytical solution coefficients
    L1pdd = L1p+k3+k4; L1pd = L1pdd + k56; L1mdd = L1m+k3+k4; L1md = L1mdd + k56;
    L2pdd = L2p+k3; L2pd = L2pdd + k56; L2mdd = L2m+k3; L2md = L2mdd + k56;
    L3pdd = L3p+k4; L3pd = L3pdd + k56; L3mdd = L3m+k4; L3md = L3mdd + k56;

    c1p = k5*L1pd - k6*L1mdd; c1m = k6*L1pdd - k5*L1md;
    c2p = k5*L2pd - k6*L2mdd; c2m = k6*L2pdd - k5*L2md;
    c3p = k5*L3pd - k6*L3mdd; c3m = k6*L3pdd - k5*L3md;

    % Terms associated with L1, L2, and L3
    T1 = (c1p*exp(L1p*t) + c1m*exp(L1m*t)) / (DL1*k56);
    T2 = (c2p*exp(L2p*t) + c2m*exp(L2m*t)) / (DL2*k56);
    T3 = (c3p*exp(L3p*t) + c3m*exp(L3m*t)) / (DL3*k56);

    % Analytical solutions for single- and double-labeled yields
    Ys = 1 - (1/r)*T1 - ((r-1)/(2*r))*T2 - ((r-1)/(2*r))*T3;
    Yd = 1 + T1 - T2 - T3;
    
    % Set output: Single: Ys; Double: Yd; Both: [Ys, Yd]
    switch mode
        case 'Double'
            Ys = Yd;
    end
end