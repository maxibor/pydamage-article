import streamlit as st
import pydamage
import statsmodels
import pkg_resources
import pandas as pd
import numpy as np
import pickle
import gzip
import pathlib

@st.cache()
def load_model():
    """Returns the gml model"""
    this_dir = str(pathlib.Path(__file__).parent.absolute())
    model_path = this_dir+"/../models/accuracy_model_v2_python.pickle.gz"
    with gzip.open(model_path, "rb") as mod:
        return pickle.load(mod)


def prepare_df(coverage, contiglength, damage):
    """Prepare dataframe from input parameters

    Args:
        coverage (float): Mean coverage of contig
        contiglength (int): Length of contig
        damage (float): Damage on 5' end of contig

    Returns:
        [pandas DataFrame]: parameters as df
    """
    var_dict = {
        "coverage": [float(coverage)],
        "contiglength": [float(contiglength)],
        "damage": [float(damage)],
    }
    pd_df = pd.DataFrame(var_dict)
    return pd_df


def fit_model(df, model):
    """Fit GLM model to data
    Args:
        df (pandas DataFrame): prepared pydamage results
        model (pypmml model): GLM accuracy model
    """
    return model.predict(df)[0]


def interactive_predict(model):
    intro = """
# [PyDamage](https://github.com/maxibor/pydamage) Accuracy prediction 

Interactive exploration of the effect of variables on PyDamage Predicted Accuracy.

<img src="https://raw.githubusercontent.com/maxibor/pydamage-article/master/plots/Predicted_Accuracy.png" alt="Simulation scheme" width="100%">  

**Figure 4**: Predicted model accuracy of simulated data.  The grey title box above each panel is the simulated damage frequency on the 5' end. Light blue indicates improved model accuracy, with parameter combinations resulting in better than 50% accuracy are outlined in green. 

See [preprint](https://www.biorxiv.org/content/10.1101/2021.03.24.436838v1) for details
    """
    st.markdown(intro, unsafe_allow_html=True)

    coverage = st.number_input(label='Contig Mean Coverage', min_value=0.0, max_value=99999999999.0, value=25.3, step=0.01)
    contiglength = st.number_input(label='Contig Length (bp)', min_value=0, max_value=99999999999, value=10000, step=1)
    damage = st.number_input(label="Damage on 5' end", min_value=0.0, max_value=0.99, value=0.3, step=0.01)

    df = prepare_df(coverage, contiglength, damage)

    pred_acc = fit_model(df, model)

    st.markdown(f"**Predicted accuracy:** `{round(pred_acc,2)}`")

if __name__ == '__main__':
    model = load_model()
    interactive_predict(model)