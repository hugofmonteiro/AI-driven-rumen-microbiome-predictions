import pandas as pd
import numpy as np
import sklearn
import statsmodels.api as sm
import sklearn.metrics as metrics
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, explained_variance_score
from sklearn.preprocessing import StandardScaler
import sklearn.feature_selection as skfs
from sklearn.metrics import accuracy_score, r2_score
from mlxtend.feature_selection import SequentialFeatureSelector
from tabulate import tabulate
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

# Importing the training and test datasets for debuging
microbiomeRFI_Variables = pd.read_csv('ALL_complete_without_repeatedAnimals.csv')
microbiomeRFI_Variables.columns.to_list()

# Dataset containing microbiome variables at the Phylum, Class, Order, Family and Genus level
# Data contain variables as microbial counts, relative abundance, and centered-log ratio normalization structure
# Remove from drop Parity, ExpN, NESec, MBW, and BEC when predicting DMI, MFE, and MPE
X = microbiomeRFI_Variables.drop(['Obs', 'ID', 'FDAT', 'Site', 'Exp', 'ExpN', 'Lact', 'Parity', 'ParityN', 'FirstDIM', 'LastDIM', 'DOBS', 'GFE', 'DMI', 'DMI_tercile', 'RFI', 'RMFE','RMPE','RMLE','PCoA1', 'PCoA2', 'DietNEL', 'NELI', 'Milk', 'Milk_tercile', 'Milk305', 'PCTF', 'MFPE', 'PCTF_MBW', 'PCTF_BEC', 'PCTP', 'PCTL', 'MFE', 'MPE', 'MLE', 'FatY', 'FatY_tercile', 'FatY_quartile', 'ProtY', 'LactY', 'ECM', 'NESec', 'BW', 'BWC', 'BCS', 'MBW' ,'BEC', 'gMilk', 'gFat', 'gProt', 'PL', 'SCS', 'DPR', 'CCR', 'HCR', 'LIV', 'GL', 'gRFI', 'MFV', 'DAB', 'KET', 'MAS', 'MET', 'RPL', 'EFC', 'HLV', 'FS', 'NM$', 'FM$', 'CM$', 'GM$', 'rMilk', 'rFat', 'rProt', 'rPL', 'rSCS', 'rDPR', 'CrCR', 'rHCR', 'rLIV', 'rGL', 'rRFI', 'rMFV', 'rDAB', 'rKET', 'rMAS', 'rMET', 'rRPL', 'rEFC', 'rHLV', 'rFS', 'rNM$', 'rFM$', 'rCM$', 'rGM$', 'RFI_half', 'GFE_groups', 'DMI_groups', 'RFI_0L1Me2Mo_groups', 'RFI_quartile',  'RFI_percentile', 'DietNEL_groups', 'NELI_groups', 'Milk_groups', 'PCTF_groups', 'FatY_PCTF_groups', 'FatY_DMI_groups', 'PCTP_groups', 'ProtY_PCTP_groups', 'ProtY_DMI_groups', 'PCTL_groups', 'LactY_PCTL_groups', 'LactY_DMI_groups', 'MFE_groups', 'MPE_groups', 'MLE_groups', 'FatY_groups', 'ProtY_groups', 'LactY_groups', 'ECM_groups', 'NESec_groups', 'BW_groups', 'MBW_groups', 'BWC_groups', 'BCS_groups', 'BEC_groups', 'gMilk_groups', 'gFat_groups', 'gProt_groups', 'PL_groups', 'SCS_groups', 'DPR_groups', 'CCR_groups', 'HCR_groups', 'LIV_groups', 'GL_groups', 'RFI_groups', 'MFV_groups', 'DAB_groups', 'KET_groups', 'MAS_groups', 'MET_groups', 'RPL_groups', 'EFC_groups', 'HLV_groups', 'FS_groups', 'NM$_groups', 'FM$_groups', 'CM$_groups', 'GM$_groups', 'rMilk_groups', 'rFat_groups', 'rProt_groups', 'rPL_groups', 'rSCS_groups', 'rDPR_groups', 'CrCR_groups', 'rHCR_groups', 'rLIV_groups', 'rGL_groups', 'rRFI_groups', 'rMFV_groups', 'rDAB_groups', 'rKET_groups', 'rMAS_groups', 'rMET_groups', 'rRPL_groups', 'rEFC_groups', 'rHLV_groups', 'rFS_groups', 'rNM$_groups', 'rFM$_groups', 'rCM$_groups', 'rGM$_groups'], axis=1)

# Correct when predicting MFE and MPE
X['Parity'] = X['Parity'].replace({'Mult': 1, 'Prim': 0})

y = microbiomeRFI_Variables.loc[:,['RFI']].values

#### Standardization

import pandas as pd
import numpy as np
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, r2_score
from mlxtend.feature_selection import SequentialFeatureSelector
from sklearn.linear_model import Ridge

#For Feature I will use StandardScaler
from sklearn.preprocessing import StandardScaler
scFeatures = StandardScaler()
X=scFeatures.fit_transform(X)

# Filtering microbial variables present in just too few cows was also tested
# Filtering 0 variables different than the average of that column
X = pd.DataFrame(X)
def delete_variables(X):
    for col in X.columns:
        if X[col].dtype != object:
            col_avg = X[col].mean()
            unique_values = X[col].nunique()
            if unique_values < 0 or unique_values < abs(X[col] - col_avg).nunique():
                X = X.drop(columns=[col])
    return X

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error, r2_score
from joblib import Parallel, delayed

def calculate_mse_for_candidate_feature(candidate_feature, selected_features):
    candidate_features = selected_features.union({candidate_feature})
    candidate_X = X.iloc[:, list(candidate_features)]

    mse_scores_fold = []
    y_val_fold = np.array([])
    y_pred_fold = np.array([])
    val_indices = []

    for train_index, val_index in kfold.split(candidate_X):
        X_train, X_val = candidate_X.iloc[train_index], candidate_X.iloc[val_index]
        y_train, y_val = y[train_index], y[val_index]

        model.fit(X_train, y_train)
        y_pred = model.predict(X_val)
        mse = mean_squared_error(y_val, y_pred)

        mse_scores_fold.append(mse)
        y_val_fold = np.concatenate((y_val_fold, y_val))
        y_pred_fold = np.concatenate((y_pred_fold, y_pred))
        val_indices.extend(val_index)

    avg_mse = np.mean(mse_scores_fold)
    return candidate_feature, avg_mse, y_val_fold, y_pred_fold, val_indices

y = y.flatten()

# Define the Ridge regression model
model = Ridge(alpha=1.0)

# Define the k-fold cross-validation
kfold = KFold(n_splits=10, shuffle=True, random_state=42)

# Initialize the feature set and scores
selected_features = set()
mse_scores = []
all_y_val = np.array([])
all_y_pred = np.array([])

# Initialize the best model and best features
best_model = None
best_features = None
lowest_mse = float('inf')

# Initialize the validation indices for the best model
best_val_indices = []

# Perform the sequential feature selection with k-fold cross-validation
for feature in range(X.shape[1]):
    # Calculate the average MSE for each candidate feature in parallel
    results = Parallel(n_jobs=-1)(
        delayed(calculate_mse_for_candidate_feature)(candidate_feature, selected_features)
        for candidate_feature in range(X.shape[1])
        if candidate_feature not in selected_features
    )

    # Find the best candidate feature based on the average MSE
    best_feature, best_mse, best_y_val, best_y_pred, best_val_indices = min(results, key=lambda x: x[1])

    selected_features.add(best_feature)
    mse_scores.append(best_mse)

    # Print the best feature and its corresponding MSE for this iteration
    print(f"Selecting feature {best_feature} with MSE: {best_mse:.4f}")

    # Update the best model and best features if the current model has a lower MSE
    current_mse = mean_squared_error(best_y_val, best_y_pred)
    if current_mse < lowest_mse:
        lowest_mse = current_mse
        best_model = model
        best_features = selected_features.copy()
        best_val_indices = best_val_indices

# Create a dataset with only the validation samples for the best model
X_best_val = X.iloc[best_val_indices, list(best_features)]
y_best_val = y[best_val_indices]

# Initialize arrays to store the predictions and true values
y_pred_best = np.array([])
y_true_best = np.array([])

# Perform 10-fold cross-validation for the best model and the best feature set
for train_index, val_index in kfold.split(X_best_val):
    X_train, X_val = X_best_val.iloc[train_index], X_best_val.iloc[val_index]
    y_train, y_val = y_best_val[train_index], y_best_val[val_index]

    best_model.fit(X_train, y_train)
    y_pred = best_model.predict(X_val)

    y_pred_best = np.concatenate((y_pred_best, y_pred))
    y_true_best = np.concatenate((y_true_best, y_val))

# Calculate R² for the best model
r2_best = r2_score(y_true_best, y_pred_best)

# Print R²
print(f"R² for the best model: {r2_best:.4f}")

# Calculate MSE and RMSE for the best model
mse_best = mean_squared_error(y_true_best, y_pred_best)
rmse_best = np.sqrt(mse_best)
print(f"MSE for the best model: {mse_best:.4f}")
print(f"RMSE for the best model: {rmse_best:.4f}")

# Display the selected features and their corresponding column names
print("Selected feature for the best model:")
for feature in best_features:
    print(f"Column index: {feature}, Column name: {X.columns[feature]}")

# Set seaborn style
sns.set_style("whitegrid")

# Changing the color palette
ucd_blue = '#002855'
ucd_gold = '#FFC72C'

# Plot observed vs. predicted values for the best model
plt.figure(figsize=(10, 10))
plt.scatter(y_true_best, y_pred_best, alpha=0.7, color=ucd_blue, edgecolors='k', linewidths=0.5, s=60)
plt.xlabel('Observed RFI, kg/d', fontsize=18)
plt.ylabel('pMicrobiome RFI, kg/d', fontsize=18)

# Add a diagonal line
min_val, max_val = min(y_true_best.min(), y_pred_best.min()), max(y_true_best.max(), y_pred_best.max())
plt.plot([min_val, max_val], [min_val, max_val], color=ucd_gold, linestyle='--', lw=2)

# Add R², MSE, and RMSE to the plot
text = f'R²: {r2_best:.4f}\nMSE: {mse_best:.4f}\nRMSE: {rmse_best:.4f}'
plt.text(0.05, 0.95, text, transform=plt.gca().transAxes, fontsize=18, verticalalignment='top')

# Change plot frame colors
ax = plt.gca()
ax.spines['bottom'].set_color('k')
ax.spines['top'].set_color('k')
ax.spines['right'].set_color('k')
ax.spines['left'].set_color('k')

# Change axis label colors
ax.xaxis.label.set_color('k')
ax.yaxis.label.set_color('k')
ax.tick_params(colors='k')

# Increase the font size of the x and y tick labels
ax.tick_params(axis='x', labelsize=18)
ax.tick_params(axis='y', labelsize=18)

plt.savefig('RFI_pMicrobiome_observed_vs_predicted_Ridge.png', dpi=300, bbox_inches='tight')
plt.show()

# Create a DataFrame with the observed and predicted values
df_RFI_pred = pd.DataFrame({'Observed RFI, kg/d': y_true_best, 'Predicted RFI, kg/d': y_pred_best})

# Save the DataFrame as a CSV file
df_RFI_pred.to_csv('RFI_pMicrobiome_observed_vs_predicted.csv', index=False)

df_RFI_pred

# Create a new column with the categories
df_RFI_pred['Category'] = ['true positive' if (obs > 0) & (pred > 0) else
                           'true negative' if (obs < 0) & (pred < 0) else
                           'false positive' if (obs < 0) & (pred > 0) else
                           'false negative' for obs, pred in zip(df_RFI_pred['Observed RFI, kg/d'], df_RFI_pred['Predicted RFI, kg/d'])]

# Print the updated dataframe
print(df_RFI_pred)

# Save the DataFrame as a CSV file
df_RFI_pred.to_csv('RFI_pMicrobiome_predicted_categories.csv', index=False)

RFI_Rates = df_RFI_pred['Category'].value_counts(normalize=True) * 100

RFI_Rates

# Save the DataFrame as a CSV file
RFI_Rates.to_csv('RFI_pMicrobiome_Categories_Rates.csv', index=True)
