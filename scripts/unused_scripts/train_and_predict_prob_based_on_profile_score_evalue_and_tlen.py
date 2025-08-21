import os
import glob
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier

script_dir = os.path.dirname(os.path.realpath(__file__))
base_dir = os.path.join(script_dir, "..")

batch_num4train = 2 # It shows the number of batches whose data will be used for training a model.

def prepare_raw_df_data(df):
    """Adds a column for the target length and also calculates the log of the evalue"""
    df["tlen"] = df["target"].str.split("-", expand=True)[2].astype(int) - df["target"].str.split("-", expand=True)[1].astype(int) + 1
    df["evalue"] = df["evalue"].replace(0, 1e-100)
    df["log_evalue"] = np.log10(df["evalue"])
    return df

def score_ali_df(ml_model, ali_path, cols2read, feature_cols):
    """
    It takes a model as input along with the path to the alignment output to score.
    It also takes the important columns needed for reading the scoring dataframes, and the
    feature columns needed for scoring each hit. In the end, it will store the reranked hits in a file.
    """
    ali_df = pd.read_csv(ali_path, sep="\t", usecols=cols2read)
    ali_df = prepare_raw_df_data(ali_df)
    ali_df["ml_score"] = ml_model.predict_proba(ali_df[feature_cols])[:,1]
    ali_df.to_csv(ali_path.replace("_bitscore.tsv", f"_ml{batch_num4train}.tsv"), sep="\t", index=None)  ######### Added ml2.tsv to train a model based on two batches rather than one

def train_and_predict(tool):
    """The tool could be fs or reseek. It will train a model on the data of the first batch to predict the probability of having
    a True Positive hit based on the e-value, the target length, and the bitscore of the hit."""

    if tool == "fs":
        training_data_paths = [f"{base_dir}/tmp/alis/split_pf_seq/fs_pref_B{x}_bitscore.tsv" for x in range(1, 1 + batch_num4train)]  # The first batch is used for training the data
        to_score_paths = [f"{base_dir}/tmp/alis/split_pf_seq/fs_pref_B{x}_bitscore.tsv" for x in range(1 + batch_num4train, 17)] # Other batches are used for evaluating the new score based on tlen, profile bitscore, and the evalue of the hit
    elif tool == "reseek":
        training_data_paths = [f"{base_dir}/tmp/alis/split_pf_seq/reseek_sens_B{x}_bitscore.tsv" for x in range(1, 1 + batch_num4train)]  # The first batch is used for training the data
        to_score_paths = [f"{base_dir}/tmp/alis/split_pf_seq/reseek_sens_B{x}_bitscore.tsv" for x in range(1 + batch_num4train, 17)] # Other batches are used for evaluating the new score based on tlen, profile bitscore, and the evalue of the hit

    cols2read_test = ["query", "target", "evalue", "bitscore_rep"]  # The columns to read for scoring the hits

    features = ["bitscore_rep", "log_evalue", "tlen"]
    target = "pfam_label"

    df = pd.concat([pd.read_csv(training_data_path, sep="\t", usecols=cols2read_test) for training_data_path in training_data_paths] )
    df["pfam_label"] = (df["query"].str.split("-", expand=True)[3] == df["target"].str.split("-", expand=True)[3])  # Adds label to the training set
    df = prepare_raw_df_data(df)

    X = df[features]
    y = df[target]

    # Initialize the Random Forest Classifier
    rf_model = RandomForestClassifier(n_estimators=20, random_state=42, class_weight='balanced') #  The class_weight = balanced is used to deal with the unbalanced dataset
    # Fit the model on the training data
    rf_model.fit(X, y)
    del X, y, df  # To just release some RAM for the rest of analysis

    for path in to_score_paths:
        print(f"scoring {path}")
        score_ali_df(rf_model, path, cols2read_test, features)


if __name__ == "__main__":
    train_and_predict("fs")
    train_and_predict("reseek")
