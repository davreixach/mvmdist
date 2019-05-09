function obj = fitmvmdist(angles, nComponents, varargin)
% FITMVMDIST This function estimates the parameters of a mixture of von
%   Mises distributions using an Expectation-Maximization scheme. It uses
%   the corresponding functions provided by the VonMisesMixture class.
%   Please refer to the function headers in the class source file for
%   further information.
%
% REQUIRED INPUTS:
%   angles - Nx1 vector, containing N angular values, ranged between -pi
%       and pi.
%   nComponents - Number of mixture components that should be estimated.
%
% PARAMETERS:
%   ['MaxIter', maxIter] - Maximum number of iterations the EM-algorithm
%       should run (default = 100).
%   ['ErrorThreshold', errorThreshold] - Minimum error to  be used as a
%       stopping-criterion for the EM-algorithm during convergence testing
%       (default = 1E-4).
%   ['Replicates', replicates] - Number of replications of the parameter
%       estimation procedure. If the number of replications is greater than
%       one, the parameters of the replicate that yielded the highest
%       log-likelihood will be returned (default = 1).
%
% DEPENDS ON:
%   VonMisesMixture.m
%
% AUTHOR:
%   Copyright (c) 2016      Christopher Schymura
%                           Cognitive Signal Processing Group
%                           Ruhr-Universitaet Bochum
%                           Universitaetsstr. 150
%                           44801 Bochum, Germany
%                           E-Mail: christopher.schymura@rub.de

% Check inputs
p = inputParser();
defaultMaxIter = 100;
defaultErrorThreshold = 1E-4;
defaultReplicates = 1;
defaultWeightsMode = 'auto';

p.addRequired('Angles', ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'real', 'vector', '>=', -pi, '<=', pi}) ...
  );

p.addRequired('NComponents', ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'integer', 'scalar', 'positive'}) ...
  );

p.addParameter('MaxIter', ...
  defaultMaxIter, ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'integer', 'scalar', 'positive'}) ...
  );

p.addParameter('ErrorThreshold', ...
  defaultErrorThreshold, ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'real', 'scalar', 'nonnegative'}) ...
  );

p.addParameter('Replicates', ...
  defaultReplicates, ...
  @(x) validateattributes(x, ...
  {'numeric'}, ...
  {'integer', 'scalar', 'nonnegative'}) ...
  );

p.addParameter('WeightsMode', ...
  defaultWeightsMode, ...
  @(x) validateattributes(x, ...
  {'char'}, ...
  {'nonempty'}) ...
  );

p.parse(angles, nComponents, varargin{:});

% Initialize von Mises mixture models for all replicates
models = cell(p.Results.Replicates, 1);

% Initialize indicator variables for maximum log likelihood tracking
maxLogLik = -realmax;
bestIdx = 1;

% Perform parameter estimation
for rIdx = 1 : p.Results.Replicates
  % Initialize mixture of von Mises distribution
  model = VonMisesMixture();
  
    % Run EM and generate model
      models{rIdx} = model.fit(p.Results.Angles, p.Results.NComponents, ...
        'MaxIter', p.Results.MaxIter, ...
        'ErrorThreshold', p.Results.ErrorThreshold,...
        'WeightsMode',p.Results.WeightsMode);
  
%   switch p.Results.WeightsMode
%       case 'auto'
%           % Run EM and generate model
%           models{rIdx} = model.fit(p.Results.Angles, p.Results.NComponents, ...
%             'MaxIter', p.Results.MaxIter, ...
%             'ErrorThreshold', p.Results.ErrorThreshold);
%       case 'fixed'
%           % Run EM and generate model with equal weights
%           models{rIdx} = model.fit(p.Results.Angles, p.Results.NComponents, ...
%             'MaxIter', p.Results.MaxIter, ...
%             'ErrorThreshold', p.Results.ErrorThreshold,...
%             'Weight',ones(size(p.Results.NComponents)));
%       otherwise
%           error('WeightsMode can be either ''auto'' or ''fixed''.')
%   end

  
  % Update maximum log-likelihood
  if models{rIdx}.logLikelihood > maxLogLik
    maxLogLik = models{rIdx}.logLikelihood;
    bestIdx = rIdx;
  end
end

% Select and return best-matching model
obj = models{bestIdx};

end
