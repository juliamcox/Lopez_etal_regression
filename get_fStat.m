
function [Fstat, Fstat_noRefit] = get_fStat(X_full,X_reduced,Y,b,b_reduced,omitPreds)

yhat = cat(2,ones(size(X_full,1),1),X_full)*b; % estimated data 
sse_full = sum((Y - yhat).^2,1,'omitnan'); % sum squared error of full model
df_full = size(X_full,1)-numel(b); % degrees of freedom of full model

yhat_reduced = cat(2,ones(size(X_reduced,1),1),X_reduced)*b_reduced;
sse_reduced = sum((Y-yhat_reduced).^2,1,'omitnan'); % sum squared error of reduced model
df_reduced  = size(X_reduced,1)-numel(b_reduced); % degrees of freedom of full model

Fstat = ((sse_reduced-sse_full)/(df_reduced-df_full))/(sse_full/df_full);


% Calculate f-statistic without refitting

b_reduced = b;
b_reduced(omitPreds+1) = [];
yhat_reduced = cat(2,ones(size(X_reduced,1),1),X_reduced)*b_reduced;
sse_reduced = sum((Y-yhat_reduced).^2,1,'omitnan'); % sum squared error of reduced model
df_reduced  = size(X_reduced,1)-numel(b_reduced); % degrees of freedom of full model

Fstat_noRefit = ((sse_reduced-sse_full)/(df_reduced-df_full))/(sse_full/df_full);


