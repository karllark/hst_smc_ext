import numpy as np
from astropy.table import QTable

itab = QTable.read("tables/gor24_smc_forecor_ensemble_params.dat", format="ascii.ipac")
npts = len(itab)

# X is your data table, where the features (bump strength, C2, B3, etc) are columns and 
# individual sightlines are rows
X = np.zeros((npts, 5))
#X[:, 0] = itab["C1"]
X[:, 0] = itab["C2"]
X[:, 1] = itab["B3"]
X[:, 2] = itab["C4"]
X[:, 3] = itab["RV"].data
X[:, 4] = itab["NHI"].data / itab["EBV"].data

N_sightlines, N_features = X.shape

# normalize your data, i.e., give all features zero mean and unit variance
X_norm = (X - X.mean(axis=0)) / X.std(axis=0)

# compute covariance matrix between features; note that "@" is matrix multiply
# cov_X.shape should equal (N_features, N_features)
cov_X = np.cov(X_norm.transpose())

# decompose covariance matrix into eigenvalues/eigenvectors
# eigenvalues_X.shape should equal (N_features,)
# eigenvectors_X.shape should equal (N_features, N_features)
eigenvalues_X, eigenvectors_X = np.linalg.eig(cov_X)

# sort them in order from highest to lowest eigenvalues; now your highest eigenvalue  
# matches the eigenvector that explains most of the variance, aka your first principal component
# second highest eigenvalue is your second principal component, and etc
order = np.argsort(np.abs(eigenvalues_X))[::-1]

# reorder the eigenvalues and eigenvectors
eigenvalues_X = eigenvalues_X[order]
principal_components = eigenvectors_X[:, order]

# explained variance is simply the eigenvalue divided by sum of all eigenvalues
explained_variance = eigenvalues_X / eigenvalues_X.sum()

print(explained_variance)

for k in range(len(explained_variance)):
    print(f"Component #{k+1} explains {np.round(100*explained_variance[k], 3)} % of the variance, and is in the direction", principal_components[:, k])
