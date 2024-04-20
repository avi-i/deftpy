from sklearn.linear_model import HuberRegressor
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from adjustText import adjust_text
from sklearn.model_selection import KFold, train_test_split, cross_val_score
from sklearn.metrics import mean_absolute_error, mean_squared_error
# import seaborn

def main():

    # read in df
    df_cfm = pd.read_csv("wexler_nw_match.csv")
    # df_cfm = df_cfm.drop(df_cfm[df_cfm['Ev'] < 0].index)
    df_cfm = df_cfm[["formula", "Ev", "Eg", "vacancy", "space_group", "Eb", "Vr_max", "PT_oxide"]].reset_index(drop=True)
    df_cfm = df_cfm.dropna()
    df_cfm = df_cfm.reset_index(drop=True)
    # band gaps (needs to be done using wexler_ww.csv)
    # calculated = df_cfm["Eg"]
    # mp = df_cfm["mp_Eg"]
    # unique_formulas = df_cfm['formula'].unique()
    # colors = plt.cm.Paired(np.linspace(0, 1, len(unique_formulas)))
    # formula_color_map = {formula: color for formula, color in zip(unique_formulas, colors)}
    # # Plot scatter plot with colors mapped to formulas
    # plt.scatter(mp, calculated)
    # plt.plot(mp, mp, "k--")
    # plt.gca().set_aspect('equal', adjustable='box')
    # texts = []
    # for x, y, s in zip(mp, calculated, df_cfm["space_group"]):
    #     texts.append(plt.text(x, y, s, size=10))
    # adjust_text(texts, arrowprops=dict(arrowstyle="-", color="k", lw=0.5))
    # plt.scatter(df_cfm["mp_Eg"], df_cfm["Eg"], c=df_cfm["formula"].map(formula_color_map))
    # plt.title("band gaps from MP vs. DFT")
    # plt.xlabel("materials project")
    # plt.ylabel("DFT calculated")
    # # plt.legend(formula_color_map, loc=(1.1, 0.5))
    # handles = []
    # labels = []
    # for key, value in formula_color_map.items():
    #     handles.append(plt.Line2D([], [], marker='o', color=value, linestyle='None', markersize=10))
    #     labels.append(key)
    # plt.legend(handles, labels, loc=(1.05, 0.1))
    # plt.savefig("wexler_mp_calc_Eg.png")
    # plt.show()

    # create definitions
    cfm = HuberRegressor()
    X = df_cfm[["Vr_max"]]
    X = df_cfm[["Eb", "Vr_max"]]
    X = df_cfm[["Eb", "Vr_max", "Eg"]]
    # X = df_cfm[["Eb_sum", "Vr_max", "mp_Eg"]]
    # X = df_cfm[["Eb_sum", "Vr_max", "mp_Eg", "mp_ehull"]]
    y = df_cfm["Ev"]
    cfm.fit(X, y)
    y_pred = cfm.predict(X)
    cfm.score(X, y)

    kf = KFold(n_splits=5, shuffle=True, random_state=21)
    scores = cross_val_score(cfm, X, y, scoring='neg_mean_absolute_error', cv=kf)
    coefs = cfm.coef_
    for score in scores:
        print('score for this fold is ', score)
        print('coefficients', cfm.coef_)

    # plt.hist(scores, bins=20)
    # plt.show()

    mean = np.mean(scores)
    std = np.std(scores)
    rmse = mean_squared_error(y, y_pred, squared=False)

    unique_formulas = df_cfm['formula'].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_formulas)))
    formula_color_map = {formula: color for formula, color in zip(unique_formulas, colors)}
    unique_periods = df_cfm["PT_oxide"].unique()
    marker = {f"{unique_periods[0]}": 'o', f"{unique_periods[1]}": '^', f"{unique_periods[2]}": 's'}

    plt.plot([1, 9], [1, 9], "k--")

    for _, row in df_cfm.iterrows():
        oxide = row["PT_oxide"]
        plt.scatter(y_pred[_], y[_], color=(formula_color_map[row["formula"]]), marker=marker[oxide])

    # plt.scatter(y_pred, y)
    plt.xlim(min(y_pred) - 1, max(y_pred) + 1)
    plt.ylim(min(y) - 1, max(y) + 1)
    plt.gca().set_aspect('equal', adjustable='box')
    equation = "$E_v = {:.2f} {:+.2f} E_b {:+.2f} V_r {:+.2f} E_g$".format(cfm.intercept_, cfm.coef_[0],
                                                                               cfm.coef_[1], cfm.coef_[2])
    mae = np.mean(np.abs(y - y_pred))
    # oxides = f"$MO_x$, $ABO_x$"
    oxides = f"$MO_x$"
    plt.text(5, 0, f"MAE = {mae:.2f}", fontsize=10)
    # plt.text(7, 4, oxides, fontsize=12)
    plt.text(5, -0.5, f"n = {len(y)}", fontsize=10)
    plt.text(5, -1.0, f"mean = {-1 * mean:.2f}", fontsize=10)
    plt.text(5, -1.5, f"std = {std:.2f}", fontsize=10)
    plt.xlabel(str(equation))
    texts = []
    # for x, y, s in zip(y_pred, y, df_cfm["space_group"]):
    #         texts.append(plt.text(x, y, s, size=6))
    # adjust_text(texts, arrowprops=dict(arrowstyle="-", color="k", lw=0.5))
    plt.ylabel(f"Wexler $E_v$")
    plt.title("CFM for binaries with non-weighted CN")
    handles = []
    labels = []
    for key, value in formula_color_map.items():
        handles.append(plt.Line2D([], [], marker='o', color=value, linestyle='None', markersize=8))
        labels.append(key)
    # plt.legend(handles, labels, loc=(1.05, 0.1))
    for key, value in marker.items():
        handles.append(plt.Line2D([], [], marker=value, linestyle='None', markersize=8))
    for x in ['alkaline earth', 'alkali', 'group 13']:
        labels.append(x)
    plt.legend(handles, labels, loc=(1.05, 0.08))
    plt.savefig("wexler_Eg_frac_outlier_nw.pdf")
    plt.show()

if __name__ == "__main__":
    main()