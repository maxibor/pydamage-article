import streamlit as st
import pkg_resources
import os
from pypmml import Model
import pickle

def prepare_dict(coverage, contiglength, damage):
    """Prepare dataframe from input parameters

    Args:
        coverage (float): Mean coverage of contig
        contiglength (int): Length of contig
        damage (float): Damage on 5' end of contig

    Returns:
        [pandas DataFrame]: parameters as df
    """
    var_dict = {
        "actual_cov": float(coverage),
        "damage": float(damage),
        "contiglength": int(contiglength)
    }
    print(var_dict)
    return var_dict


def predict_model(model, var_dict):
    """Fit GLM model to data
    Args:
        var_dict (var_dict): prepared pydamage results
    """
    print(model.predict({'actual_cov':44.0, 'damage':0.01, 'contiglength':800}))
    return model.predict(var_dict)['Predicted_sig']


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

    var_dict = prepare_dict(coverage, contiglength, damage)

    pred_acc = predict_model(model, var_dict)

    st.markdown(f"**Predicted accuracy:** `{round(pred_acc,2)}`")

def test(model):

    print("test", model.predict({'actual_cov':44.0, 'damage':0.01, 'contiglength':800}))

if __name__ == '__main__':
    dirname = os.path.dirname(__file__)
    model_path = dirname + "/../models/pydamage_glm_model.pmml"
    model = Model.load(model_path)
    print("test2", model.predict({'actual_cov':23.0, 'damage':0.3, 'contiglength':45678}))
    print(model.inputNames)
    test(model)
    interactive_predict(model)
