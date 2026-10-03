"""
Statistical analysis module for SQUID coercivity and remanence datasets.

Implements OLS joint linear regression, two-way categorical ANOVA,
confidence interval calculation, and summary statistics using numpy and scipy.stats.
"""

import numpy as np
import scipy.stats


def run_ols_joint_regression(diameters, coercivities, batch_indicators):
    """
    Fit joint linear regression across batches:
    H_c = beta_0 + beta_batch1 * I_batch1 + beta_slope * (1/R)
    where R = D / 2 in micrometers (so 1/R = 2/D in um^-1).
    
    Parameters
    ----------
    diameters : list or np.ndarray
        Sphere diameters in micrometers (e.g. [3, 5, 8, 10, ...]).
    coercivities : list or np.ndarray
        Measured H_c in Tesla (unrounded).
    batch_indicators : list or np.ndarray
        1 for Batch 1, 0 for Batch 2.
        
    Returns
    -------
    results : dict
        Detailed regression results and bounds.
    """
    y = np.array(coercivities, dtype=float)
    d = np.array(diameters, dtype=float)
    b1 = np.array([1.0 if (x == 1 or x is True or '1' in str(x)) else 0.0 for x in batch_indicators], dtype=float)
    inv_R = 2.0 / d
    const = np.ones(len(y), dtype=float)
    
    X = np.column_stack([const, b1, inv_R])
    n = len(y)
    p = X.shape[1]
    df = n - p
    
    beta, residuals, rank, s = np.linalg.lstsq(X, y, rcond=None)
    residuals = y - X @ beta
    ssr = np.sum(residuals**2)
    s2 = ssr / df
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    t_vals = beta / se
    p_vals = 2.0 * (1.0 - scipy.stats.t.cdf(np.abs(t_vals), df))
    
    # Total sum of squares for R^2
    sst = np.sum((y - np.mean(y))**2)
    r2 = 1.0 - (ssr / sst)
    
    # Delta Hc between 10 um and 3 um
    # inv_R(3 um) = 2/3, inv_R(10 um) = 2/10 = 0.2
    delta_inv_R = (2.0 / 3.0) - (2.0 / 10.0)  # 7/15 = 0.466667 um^-1
    delta_Hc = beta[2] * delta_inv_R
    se_delta = se[2] * delta_inv_R
    t_crit = scipy.stats.t.ppf(0.975, df)
    ci_lower = delta_Hc - t_crit * se_delta
    ci_upper = delta_Hc + t_crit * se_delta
    p_delta = p_vals[2]
    
    results = {
        'const': beta[0],
        'intercept': beta[0],
        'const_se': se[0],
        'const_p': p_vals[0],
        'b1_offset': beta[1],
        'batch1_offset': beta[1],
        'beta_batch1': beta[1],
        'b1_se': se[1],
        'b1_p': p_vals[1],
        'slope': beta[2],
        'slope_inv_R': beta[2],
        'slope_inv_r': beta[2],
        'slope_se': se[2],
        'slope_p': p_vals[2],
        'p_value_slope': p_vals[2],
        'r2': r2,
        'df': df,
        'residual_sd': np.sqrt(s2),
        'delta_Hc_10_to_3': delta_Hc,
        'delta_hc_10_to_3': delta_Hc,
        'delta_Hc_se': se_delta,
        'ci_95_lower': ci_lower,
        'ci_95_upper': ci_upper,
        'p_value': p_delta,
        'relative_bound_pct': (ci_upper / np.mean(y)) * 100.0
    }
    return results


def run_categorical_anova(diameters, coercivities, batch_indicators):
    """
    Fit two-way model with diameter as categorical factor + batch offset:
    Full model: H_c ~ const + I_batch1 + I_5um + I_8um + I_10um (ref: 3 um)
    Reduced model: H_c ~ const + I_batch1
    
    Returns F(3,3), p-value, and diameter effects relative to 3 um.
    """
    y = np.array(coercivities, dtype=float)
    d = np.array(diameters, dtype=float)
    b1 = np.array([1.0 if (x == 1 or x is True or '1' in str(x)) else 0.0 for x in batch_indicators], dtype=float)
    const = np.ones(len(y), dtype=float)
    
    d5 = (d == 5.0).astype(float)
    d8 = (d == 8.0).astype(float)
    d10 = (d == 10.0).astype(float)
    
    X_full = np.column_stack([const, b1, d5, d8, d10])
    df_full = len(y) - X_full.shape[1]  # 8 - 5 = 3
    beta_full = np.linalg.lstsq(X_full, y, rcond=None)[0]
    ssr_full = np.sum((y - X_full @ beta_full)**2)
    s2_full = ssr_full / df_full
    
    X_red = np.column_stack([const, b1])
    df_red = len(y) - X_red.shape[1]  # 8 - 2 = 6
    beta_red = np.linalg.lstsq(X_red, y, rcond=None)[0]
    ssr_red = np.sum((y - X_red @ beta_red)**2)
    
    f_stat = ((ssr_red - ssr_full) / (df_red - df_full)) / (ssr_full / df_full)
    p_val = 1.0 - scipy.stats.f.cdf(f_stat, df_red - df_full, df_full)
    
    return {
        'F_stat': f_stat,
        'F_statistic': f_stat,
        'df1': df_red - df_full,
        'df2': df_full,
        'df_diam': df_red - df_full,
        'df_res': df_full,
        'p_value': p_val,
        'effect_5um_vs_3um': beta_full[2],
        'effect_8um_vs_3um': beta_full[3],
        'effect_10um_vs_3um': beta_full[4],
        'residual_sd': np.sqrt(s2_full)
    }


def compute_batch_summary(values):
    """
    Compute mean and sample standard deviation (ddof=1).
    """
    arr = np.array(values, dtype=float)
    return {
        'mean': float(np.mean(arr)),
        'std': float(np.std(arr, ddof=1)),
        'min': float(np.min(arr)),
        'max': float(np.max(arr))
    }
