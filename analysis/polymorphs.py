import pandas as pd

def filter_df(df):
    # Grouping by column 'A' and checking if any group has more than one unique value in column 'B' -> Ef or Eg
    mask = df.groupby('formula')['formation_energy'].transform('nunique') > 1
    # Filtering DataFrame based on the mask
    filtered_df = df[mask]
    return filtered_df

def main():
    df = pd.read_csv('kumagai_full_Eb_BVS_Vr.csv')
    df_filtered = filter_df(df)
    print(df)
    print(df_filtered)

if __name__ == "__main__" :
    main()