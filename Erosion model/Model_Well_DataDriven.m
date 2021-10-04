close all
clear all
clc

%liquid flowrate
Ql = [2:10]'; % [L/min]

% DP [mbar]
% Use these models to calculate the max(new) / min(old) value of the
% DP in [mbar] for different flowrates
DPsampleNew = 0.2788*Ql.^2 + 1.143*Ql - 2.3831; 
DPsampleOld = 0.2693*Ql.^2 + 0.0504*Ql + 0.1256;

% Time [min]
RUL = -3*Ql + 51.333;

% model with second order terms 

% dPModel2 = a + b*tData + c*QlData + d*tData.*QlData + e*QlData.^2;

% dPModel2 = -4.45057218353704 + 0.150265436491970*tData + 1.88054326803212*QlData -0.0503621574689086*tData.*QlData + 0.213170980385915*QlData.^2;

% model with interaction terms 
%dPModel = a + b*tData + c*QlData + d*tData.*QlData;
%dPModel = -12.6939161211782 + 0.212753805754387*tData + 4.70493852121859*QlData -0.0618564737866064*tData.*QlData;
%a = dPModel2;
%% checking the behavior of the max/min DP
% the maximum dP value indicates a new probe
% the minimum dP value indicates an old probe
figure(1)
plot(Ql,DPsampleNew,'-og',Ql,DPsampleOld,'-dr')
legend({'new','old'},'Location','northwest')
xlabel('Ql [L/min]') ; ylabel('dP [mbar]') ; title('New vs. Old Sample')

%% Extrapolating
%First indicate individual columns
t = [zeros(length(Ql),1);RUL];
x = [Ql;Ql];
y = [DPsampleNew;DPsampleOld];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use the function 'griddata' for extrapolating %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first creat a grid (first coordinate is time, second is the DP)
[TI,YI] = meshgrid(0:1:46,0:1:35) ;
% extrapolate based on the individual data
XI = griddata(t,y,x,TI,YI);

% showing the extrapolated surface
figure(2)
plot3(TI,XI,YI,'Marker','o' )
xlabel('time [min]') ; ylabel('Ql [L/min]') ; zlabel('dP [mbar]') ; title('Extrapolated data (raw)')

%% data cleaning

% excluding from extrapolation: NaN values & points with erroneous behavior
% We first replace them by zeros
for ii = 1:36
    for jj = 1:47
        if isnan(XI(ii,jj))
            XI(ii,jj) = 0;
            TI(ii,jj) = 0;
            YI(ii,jj) = 0;
        end
    end
end

% Than we replace the pressure values that indicate sudden spikes (see Figure 2) 
% Cue to RUL and the extrapolation, DP values spike after the RUL is exhausted
for jj = 1:47
    for ii = 36:-1:2
        if XI(ii,jj) - XI(ii - 1,jj) < 0.01
            XI(ii,jj) = 0;
            TI(ii,jj) = 0;
            YI(ii,jj) = 0;
        end
    end
end

% size of the grid
dataLength = size(TI,1)*size(TI,2);

% reshaping data into vectors
QlData = reshape(XI,1,dataLength);
tData = reshape(TI,1,dataLength);
dPData = reshape(YI,1,dataLength);

% finding where the zeros are (from the data cleaning step)
[I,J] = find(QlData < 0.001);

% emptying the entries with zero
QlData(J) = [];
tData(J) = [];
dPData(J) = [];

% plot
figure(3)
plot3(tData,QlData,dPData,'o')
xlabel('time [min]') ; ylabel('Ql [L/min]') ; zlabel('dP [mbar]') ; title('Extrapolated data (clean)')

%% preparing data bundle
DataReg = [tData', QlData', dPData'];

%% Model with second order terms
%%%%%%%%%%%%%%%%%%%%%
% Model with Casadi %
%%%%%%%%%%%%%%%%%%%%%
import casadi.*

% Declare variables to the model (parameters and controls)
u_t = MX.sym('u_t');
u_Q = MX.sym('u_Q');
controls = [u_t;u_Q];

%N.B.: model structure was chosen based on intuition
a = MX.sym('a');
b = MX.sym('b');
c = MX.sym('c');
d = MX.sym('d');
e = MX.sym('e');
params = [a;b;c;d;e];

J = a*u_t + b*u_t*u_Q + c*u_Q + d*u_Q^2 + e;

% Create a function that simulates the quadratic surface
quadFun = Function('quad',{controls, params},{J});
%%%%%%%%%%%%%%
% Simulation %
%%%%%%%%%%%%% %
%simulating the system for all the inputs at the same time
N = size(DataReg,1);%number of points used in the estimation procedure - changes
all_samples = quadFun.map(N);

Y_symbolic = all_samples(DataReg(:,1:2)', repmat(params,1,N));
%u is tranpose: each sample corresponds to a column!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimating model parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = DataReg(:,3) - Y_symbolic';

%formalizing nlp structure and nlp problem
nlp = struct('x', params, 'f', 0.5*dot(e,e));
solver = nlpsol('solver','ipopt', nlp);

a0 = 0.149791128163722;
b0 = -0.0444958219035729;
c0 = 4.19836498078994;
d0 = 4.19836498078994;
e0 = -11.0701010584998;
xGuess = [a0;b0;c0;d0;e0];
% Solve
sol = solver('x0',xGuess,'lbx',-100*ones(5,1),'ubx',100*ones(5,1));
%sol = solver('x0',p0);
par = full(sol.x);

%%%%%%%%%%%%%%%%
% Casadi Model %
%%%%%%%%%%%%%%%%
dPModel2 = par(5) + par(1)*tData + par(3)*QlData + par(2)*tData.*QlData + par(4)*QlData.^2;

figure(4)
subplot(2,1,1)
    plot3(tData,QlData,dPData,'bo')
    hold on 
    plot3(tData,QlData,dPModel2,'rx')
    legend({'Data','Model'})
    xlabel('time [min]') ; ylabel('Ql [L/min]') ; zlabel('dP [mbar]'); title('Model w/ 2nd order terms')

subplot(2,1,2)
    plot3(tData,QlData,abs((dPData - dPModel2)./dPData),'kx')
    legend({'residual'})
    xlabel('time [min]') ; ylabel('Ql [L/min]') ; zlabel('[%]');

%% Stepwise linear model

%Without PCA
[StepwiseModel, StepwiseValidationRMSE] = trainRegressionModel(DataReg);

dPModel = StepwiseModel.predictFcn(DataReg(:,1:2));

figure(5)
subplot(2,1,1)
    plot3(tData,QlData,dPData,'bo')
    hold on 
    plot3(tData,QlData,dPModel','rx')
    legend({'Data','Model'})
    xlabel('time [min]') ; ylabel('Ql [L/min]') ; zlabel('dP [mbar]'); title('Model w/ interaction terms')

subplot(2,1,2)
    plot3(tData,QlData,abs((dPData - dPModel')./dPData),'kx')
    legend({'residual'})
    xlabel('time [min]') ; ylabel('Ql [L/min]') ; zlabel('[%]');

%% Saving Model   
name = 'ModelErosionRig';
save(name,'par','StepwiseModel','tData','QlData','dPData');


function [trainedModel, validationRMSE] = trainRegressionModel(trainingData)
% [trainedModel, validationRMSE] = trainRegressionModel(trainingData)
% returns a trained regression model and its RMSE. This code recreates the
% model trained in Regression Learner app. Use the generated code to
% automate training the same model with new data, or to learn how to
% programmatically train models.
%
%  Input:
%      trainingData: a matrix with the same number of columns and data type
%       as imported into the app.
%
%  Output:
%      trainedModel: a struct containing the trained regression model. The
%       struct contains various fields with information about the trained
%       model.
%
%      trainedModel.predictFcn: a function to make predictions on new data.
%
%      validationRMSE: a double containing the RMSE. In the app, the
%       History list displays the RMSE for each model.
%
% Use the code to train the model with new data. To retrain your model,
% call the function from the command line with your original data or new
% data as the input argument trainingData.
%
% For example, to retrain a regression model trained with the original data
% set T, enter:
%   [trainedModel, validationRMSE] = trainRegressionModel(T)
%
% To make predictions with the returned 'trainedModel' on new data T2, use
%   yfit = trainedModel.predictFcn(T2)
%
% T2 must be a matrix containing only the predictor columns used for
% training. For details, enter:
%   trainedModel.HowToPredict

% Auto-generated by MATLAB on 16-Jun-2021 18:00:44


% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3'});

predictorNames = {'column_1', 'column_2'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_3;
isCategoricalPredictor = [false, false];

% Train a regression model
% This code specifies all the model options and trains the model.
concatenatedPredictorsAndResponse = predictors;
concatenatedPredictorsAndResponse.column_3 = response;
linearModel = stepwiselm(...
    concatenatedPredictorsAndResponse, ...
    'linear', ...
    'Upper', 'interactions', ...
    'NSteps', 1000, ...
    'Verbose', 0);

% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
linearModelPredictFcn = @(x) predict(linearModel, x);
trainedModel.predictFcn = @(x) linearModelPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedModel.LinearModel = linearModel;
trainedModel.About = 'This struct is a trained model exported from Regression Learner R2019b.';
trainedModel.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 2 columns because this model was trained using 2 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appregression_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', {'column_1', 'column_2', 'column_3'});

predictorNames = {'column_1', 'column_2'};
predictors = inputTable(:, predictorNames);
response = inputTable.column_3;
isCategoricalPredictor = [false, false];

% Perform cross-validation
KFolds = 5;
cvp = cvpartition(size(response, 1), 'KFold', KFolds);
% Initialize the predictions to the proper sizes
validationPredictions = response;
for fold = 1:KFolds
    trainingPredictors = predictors(cvp.training(fold), :);
    trainingResponse = response(cvp.training(fold), :);
    foldIsCategoricalPredictor = isCategoricalPredictor;
    
    % Train a regression model
    % This code specifies all the model options and trains the model.
    concatenatedPredictorsAndResponse = trainingPredictors;
    concatenatedPredictorsAndResponse.column_3 = trainingResponse;
    linearModel = stepwiselm(...
        concatenatedPredictorsAndResponse, ...
        'linear', ...
        'Upper', 'interactions', ...
        'NSteps', 1000, ...
        'Verbose', 0);
    
    % Create the result struct with predict function
    linearModelPredictFcn = @(x) predict(linearModel, x);
    validationPredictFcn = @(x) linearModelPredictFcn(x);
    
    % Add additional fields to the result struct
    
    % Compute validation predictions
    validationPredictors = predictors(cvp.test(fold), :);
    foldPredictions = validationPredictFcn(validationPredictors);
    
    % Store predictions in the original order
    validationPredictions(cvp.test(fold), :) = foldPredictions;
end

% Compute validation RMSE
isNotMissing = ~isnan(validationPredictions) & ~isnan(response);
validationRMSE = sqrt(nansum(( validationPredictions - response ).^2) / numel(response(isNotMissing) ));
end
